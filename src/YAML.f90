!----------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parser for YAML files.
!> @details Module converts a YAML input file to an equivalent YAML flow style which is then parsed.
!----------------------------------------------------------------------------------------------------
module YAML
  use prec
  use misc
  use IO
  use types
#ifdef FYAML
  use system_routines
#endif

  implicit none(type,external)
  private

  public :: &
    YAML_init, &
    YAML_str_asList, &
    YAML_str_asDict

#ifdef FYAML
  interface

    subroutine to_flow_C(flow,length_flow,mixed) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_LONG, C_CHAR, C_PTR
      implicit none(type,external)

      type(C_PTR), intent(out) :: flow
      integer(C_LONG), intent(out) :: length_flow
      character(kind=C_CHAR), dimension(*), intent(in) :: mixed
    end subroutine to_flow_C

  end interface
#endif


contains

!--------------------------------------------------------------------------------------------------
!> @brief Do sanity checks.
!--------------------------------------------------------------------------------------------------
subroutine YAML_init()

  print'(/,1x,a)', '<<<+-  YAML init  -+>>>'
#ifdef FYAML
  print'(/,1x,a)', 'libfyaml powered'
#else
  call YAML_selfTest()
#endif

end subroutine YAML_init


!--------------------------------------------------------------------------------------------------
!> @brief Parse a YAML string with list at root into a structure of nodes.
!> @details The string needs to end with a newline (unless using libfyaml).
!--------------------------------------------------------------------------------------------------
function YAML_str_asList(str) result(list)

  character(len=*), intent(in) :: str
  type(tList), pointer :: list

  class(tNode), pointer :: node


  node => parse_flow(to_flow(str))
  list => node%asList()

end function YAML_str_asList


!--------------------------------------------------------------------------------------------------
!> @brief Parse a YAML string with dict at root into a structure of nodes.
!> @details The string needs to end with a newline (unless using libfyaml).
!--------------------------------------------------------------------------------------------------
function YAML_str_asDict(str) result(dict)

  character(len=*), intent(in) :: str
  type(tDict), pointer :: dict

  class(tNode), pointer :: node


  node => parse_flow(to_flow(str))
  dict => node%asDict()

end function YAML_str_asDict


!--------------------------------------------------------------------------------------------------
!> @brief Read a string in flow style and store it in the form of dictionaries, lists, and scalars.
!> @details A node-type pointer can either point to a dictionary, list, or scalar type entities.
!--------------------------------------------------------------------------------------------------
recursive function parse_flow(YAML_flow) result(node)

  character(len=*), intent(in) :: YAML_flow                                                         !< YAML file in flow style
  class(tNode), pointer        :: node

  class(tNode), pointer :: &
    myVal
  character(len=:), allocatable :: &
    flow_string, &
    key
  integer(pI64) :: &
    e, &                                                                                            ! end position of dictionary or list
    s, &                                                                                            ! start position of dictionary or list
    d                                                                                               ! position of key: value separator (':')


  flow_string = trim(adjustl(YAML_flow))
  if (len_trim(flow_string,pI64) == 0_pI64) then
    node => emptyDict
    return
  elseif (flow_string(1:1) == '{') then                                                             ! start of a dictionary
    e = 1_pI64
    allocate(tDict::node)
    do while (e < len_trim(flow_string,pI64))
      s = e
      d = s + scan(flow_string(s+1_pI64:),':',kind=pI64)
      e = d + find_end(flow_string(d+1_pI64:),'}')
      key = trim(adjustl(flow_string(s+1_pI64:d-1_pI64)))
      if (quotedStr(key)) key = key(2_pI64:len(key,kind=pI64)-1_pI64)
      myVal => parse_flow(flow_string(d+1_pI64:e-1_pI64))                                           ! parse items (recursively)

      select type (node)
        class is (tDict)
          call node%set(key,myVal)
      end select
    end do
  elseif (flow_string(1:1) == '[') then                                                             ! start of a list
    e = 1_pI64
    allocate(tList::node)
    do while (e < len_trim(flow_string,pI64))
      s = e
      e = s + find_end(flow_string(s+1_pI64:),']')
      myVal => parse_flow(flow_string(s+1_pI64:e-1_pI64))                                           ! parse items (recursively)

      select type (node)
        class is (tList)
          call node%append(myVal)
      end select
    end do
  else                                                                                              ! scalar value
    allocate(tScalar::node)
      select type (node)
        class is (tScalar)
          if (quotedStr(flow_string)) then
            node = trim(adjustl(flow_string(2_pI64:len(flow_string,kind=pI64)-1_pI64)))
          else
            node = trim(adjustl(flow_string))
          end if
      end select
  end if

end function parse_flow


!--------------------------------------------------------------------------------------------------
!> @brief Find location of chunk end: ',' '}', or ']'.
!> @details leaves nested lists ( '[...]' and dicts '{...}') intact
!--------------------------------------------------------------------------------------------------
integer(pI64) function find_end(str,e_char)

  character(len=*), intent(in) :: str                                                               !< chunk of YAML flow string
  character,        intent(in) :: e_char                                                            !< end of list/dict  ( '}' or ']')

  integer(pI64) :: N_sq, &                                                                          !< number of open square brackets
                   N_cu, &                                                                          !< number of open curly brackets
                                  i

  N_sq = 0_pI64
  N_cu = 0_pI64
  i = 1_pI64
  do while(i<=len_trim(str,pI64))
    if (scan(str(i:i),IO_QUOTES,kind=pI64) == 1_pI64)  i = i + scan(str(i+1:),str(i:i),kind=pI64)
    if (N_sq==0 .and. N_cu==0 .and. scan(str(i:i),e_char//',',kind=pI64) == 1_pI64) exit
    N_sq = N_sq + merge(1_pI64,0_pI64,str(i:i) == '[')
    N_cu = N_cu + merge(1_pI64,0_pI64,str(i:i) == '{')
    N_sq = N_sq - merge(1_pI64,0_pI64,str(i:i) == ']')
    N_cu = N_cu - merge(1_pI64,0_pI64,str(i:i) == '}')
    i = i + 1_pI64
  end do
  find_end = i

end function find_end


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string is enclosed with single or double quotes.
!--------------------------------------------------------------------------------------------------
logical function quotedStr(line)

  character(len=*), intent(in) :: line


  quotedStr = .false.

  if (len(line,kind=pI64) == 0) return

  if (scan(line(:1),IO_QUOTES,kind=pI64) == 1_pI64) then
    quotedStr = .true.
    if (line(len(line,kind=pI64):len(line,kind=pI64)) /= line(:1)) call IO_error(710,ext_msg=line)
  end if

end function quotedStr


#ifdef FYAML
!--------------------------------------------------------------------------------------------------
! @brief Convert all block-style YAML parts to flow style.
!--------------------------------------------------------------------------------------------------
function to_flow(mixed) result(flow)

  character(len=*), intent(in) :: mixed
  character(:,C_CHAR), allocatable :: flow

  type(C_PTR) :: str_ptr
  integer(C_LONG) :: strlen


  call to_flow_C(str_ptr,strlen,f_c_string(mixed))
  if (strlen < 1_C_LONG) call IO_error(703,ext_msg='libyfaml')
  allocate(character(len=strlen,kind=c_char) :: flow)

  block
    character(len=strlen,kind=c_char), pointer :: s
    call c_f_pointer(str_ptr,s)
    flow = s(:len(s,kind=pI64)-1_pI64)
  end block

  call free_C(str_ptr)

end function to_flow


#else
!--------------------------------------------------------------------------------------------------
! @brief Determine indentation depth.
! @details Indentation level is determined for a given block/line.
! In case of nested lists, an offset is added to determine the indent of the item block (skip
! leading dashes).
!--------------------------------------------------------------------------------------------------
integer(pI64) function indentDepth(line,offset)

  character(len=*), intent(in) :: line
  integer(pI64), optional,intent(in) :: offset


  indentDepth = verify(line,IO_WHITESPACE,kind=pI64) - 1_pI64 + misc_optional(offset,0_pI64)

end function indentDepth


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string is in flow style, i.e. starts with '{' or '['.
!--------------------------------------------------------------------------------------------------
logical function isFlow(line)

  character(len=*), intent(in) :: line


  isFlow = index(adjustl(line),'[',kind=pI64) == 1_pI64 .or. index(adjustl(line),'{',kind=pI64) == 1_pI64

end function isFlow


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string is a scalar item, i.e. starts without any special symbols.
!--------------------------------------------------------------------------------------------------
logical function isScalar(line)

  character(len=*), intent(in) :: line


  isScalar = (.not. isKeyValue(line) .and. &
              .not. isKey(line) .and. &
              .not. isListItem(line) .and. &
              .not. isFlow(line))

end function isScalar


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string is a list item, i.e. starts with '-'.
!--------------------------------------------------------------------------------------------------
logical function isListItem(line)

  character(len=*), intent(in) :: line


  isListItem = .false.
  if (len_trim(adjustl(line),pI64)> 2_pI64 .and. index(trim(adjustl(line)),'-',kind=pI64) == 1_pI64) then
    isListItem = scan(trim(adjustl(line)),' ',kind=pI64) == 2_pI64
  else
    isListItem = trim(adjustl(line)) == '-'
  end if

end function isListItem


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string contains a key-value pair of the form '<key>: <value>'.
!--------------------------------------------------------------------------------------------------
logical function isKeyValue(line)

  character(len=*), intent(in) :: line
  isKeyValue = .false.


  if ( .not. isKey(line) .and. index(clean(line),':',kind=pI64) > 0_pI64 .and. .not. isFlow(line)) then
    if (index(clean(line),': ',kind=pI64) > 0_pI64) isKeyValue = .true.
  end if

end function isKeyValue


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string contains a key without a value, i.e. it ends in ':'.
! ToDo: check whether this is safe for trailing spaces followed by a newline character
!--------------------------------------------------------------------------------------------------
logical function isKey(line)

  character(len=*), intent(in) :: line


  if (len(clean(line),kind=pI64) == 0_pI64) then
    isKey = .false.
  else
    isKey = index(clean(line),':',back=.false.,kind=pI64) == len(clean(line),kind=pI64) .and. &
            index(clean(line),':',back=.true.,kind=pI64)  == len(clean(line),kind=pI64) .and. &
            .not. isFlow(line)
  end if

end function isKey


!--------------------------------------------------------------------------------------------------
! @brief Check whether a string is a list in flow style.
!--------------------------------------------------------------------------------------------------
logical function isFlowList(line)

  character(len=*), intent(in) :: line


  isFlowList = index(adjustl(line),'[',kind=pI64) == 1_pI64

end function isFlowList


!--------------------------------------------------------------------------------------------------
! @brief Skip empty lines.
! @details Update start position in the block by skipping empty lines if present.
!--------------------------------------------------------------------------------------------------
subroutine skip_empty_lines(blck,s_blck)

  character(len=*), intent(in)     :: blck
  integer(pI64),    intent(inout)  :: s_blck

  logical :: empty


  empty = .true.
  do while (empty .and. len_trim(blck(s_blck:),pI64) /= 0_pI64)
    empty = len_trim(clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64)),pI64) == 0_pI64
    if (empty) s_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)
  end do

end subroutine skip_empty_lines


!--------------------------------------------------------------------------------------------------
! @brief Skip file header.
! @details Update start position in the block by skipping file header if present.
!--------------------------------------------------------------------------------------------------
subroutine skip_file_header(blck,s_blck)

  character(len=*), intent(in)     :: blck
  integer(pI64),    intent(inout)  :: s_blck

  character(len=:), allocatable    :: line


  line = clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64))
  if (index(adjustl(line),'%YAML',kind=pI64) == 1_pI64) then
    s_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)
    call skip_empty_lines(blck,s_blck)
    if (trim(clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64))) == '---') then
      s_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)
    else
      call IO_error(708,ext_msg = line)
    end if
  end if

end subroutine skip_file_header


!--------------------------------------------------------------------------------------------------
!> @brief Check whether a line in flow style starts and ends on the same line.
!--------------------------------------------------------------------------------------------------
logical function flow_is_closed(str,e_char)

  character(len=*), intent(in) :: str
  character,        intent(in) :: e_char                                                            !< end of list/dict  ( '}' or ']')
  integer(pI64)                :: N_sq, &                                                           !< number of open square brackets
                                  N_cu, &                                                           !< number of open curly brackets
                                  i
  character(len=:), allocatable:: line


  flow_is_closed = .false.
  N_sq = 0_pI64
  N_cu = 0_pI64
  if (e_char == ']') line = str(index(str(:),'[',kind=pI64)+1:)
  if (e_char == '}') line = str(index(str(:),'{',kind=pI64)+1:)

  do i = 1_pI64, len_trim(line,pI64)
    flow_is_closed = (N_sq==0 .and. N_cu==0 .and. scan(line(i:i),e_char,kind=pI64) == 1_pI64)
    N_sq = N_sq + merge(1_pI64,0_pI64,line(i:i) == '[')
    N_cu = N_cu + merge(1_pI64,0_pI64,line(i:i) == '{')
    N_sq = N_sq - merge(1_pI64,0_pI64,line(i:i) == ']')
    N_cu = N_cu - merge(1_pI64,0_pI64,line(i:i) == '}')
  end do

end function flow_is_closed


!--------------------------------------------------------------------------------------------------
!> @brief Return a flow-style line without line break.
!--------------------------------------------------------------------------------------------------
subroutine remove_line_break(blck,s_blck,e_char,flow_line)

  character(len=*), intent(in)               :: blck                                                !< YAML in mixed style
  integer(pI64),    intent(inout)            :: s_blck
  character,        intent(in)               :: e_char                                              !< end of list/dict  ( '}' or ']')
  character(len=:), allocatable, intent(out) :: flow_line
  logical :: line_end


  line_end = .false.
  flow_line = ''

  do while (.not. line_end)
    flow_line = flow_line//clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64))//' '
    line_end  = flow_is_closed(flow_line,e_char)
    s_blck    = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)
  end do

end subroutine remove_line_break


!--------------------------------------------------------------------------------------------------
!> @brief Return a scalar list item without line break.
!--------------------------------------------------------------------------------------------------
subroutine list_item_inline(blck,s_blck,inline,offset)

  character(len=*), intent(in)                 :: blck                                              !< YAML in mixed style
  integer(pI64),    intent(inout)              :: s_blck
  character(len=:), allocatable, intent(out)   :: inline
  integer(pI64),                 intent(inout) :: offset

  character(len=:), allocatable :: line
  integer(pI64) :: indent,indent_next


  indent = indentDepth(blck(s_blck:),offset)
  line   = clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64))
  inline = line(indent-offset+3_pI64:)
  s_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)

  indent_next = indentDepth(blck(s_blck:))

  do while (indent_next > indent)
    inline = inline//' '//trim(adjustl(clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64))))
    s_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)
    indent_next = indentDepth(blck(s_blck:))
  end do

  if (scan(inline,",",kind=pI64) > 0_pI64) inline = '"'//inline//'"'

end subroutine list_item_inline


!--------------------------------------------------------------------------------------------------
! @brief Read a line of YAML block that is already in flow style.
! @details A dict should be enclosed within '{}' for it to be consistent with the DAMASK YAML parser.
!--------------------------------------------------------------------------------------------------
recursive subroutine line_isFlow(flow,s_flow,line)

  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer(pI64),    intent(inout) :: s_flow                                                         !< start position in flow
  character(len=*), intent(in)    :: line

  integer(pI64) :: &
    s, &
    list_chunk, &
    dict_chunk


  if (index(adjustl(line),'[',kind=pI64) == 1_pI64) then
    s = index(line,'[',kind=pI64)
    flow(s_flow:s_flow) = '['
    s_flow = s_flow + 1_pI64
    do while (s < len_trim(line,pI64))
      list_chunk = s + find_end(line(s+1_pI64:),']')
      if (iskeyValue(line(s+1_pI64:list_chunk-1_pI64))) then
        flow(s_flow:s_flow) = '{'
        s_flow = s_flow + 1_pI64
        call keyValue_toFlow(flow,s_flow,line(s+1_pI64:list_chunk-1_pI64))
        flow(s_flow:s_flow) = '}'
        s_flow = s_flow + 1_pI64
      elseif (isFlow(line(s+1_pI64:list_chunk-1_pI64))) then
        call line_isFlow(flow,s_flow,line(s+1_pI64:list_chunk-1_pI64))
      else
        call line_toFlow(flow,s_flow,line(s+1_pI64:list_chunk-1_pI64))
      end if
      flow(s_flow:s_flow+1_pI64) = ', '
      s_flow = s_flow + 2_pI64
      s = s + find_end(line(s+1_pI64:),']')
    end do
    s_flow = s_flow - 1_pI64
    if (flow(s_flow-1_pI64:s_flow-1_pI64) == ',') s_flow = s_flow - 1_pI64
    flow(s_flow:s_flow) = ']'
    s_flow = s_flow + 1_pI64

  elseif (index(adjustl(line),'{',kind=pI64) == 1_pI64) then
    s = index(line,'{',kind=pI64)
    flow(s_flow:s_flow) = '{'
    s_flow = s_flow + 1_pI64
    do while (s < len_trim(line,pI64))
      dict_chunk = s + find_end(line(s+1_pI64:),'}')
      if (.not. iskeyValue(line(s+1_pI64:dict_chunk-1_pI64))) call IO_error(705,ext_msg=line)
      call keyValue_toFlow(flow,s_flow,line(s+1_pI64:dict_chunk-1_pI64))
      flow(s_flow:s_flow+1_pI64) = ', '
      s_flow = s_flow + 2_pI64
      s = s + find_end(line(s+1:),'}')
    end do
    s_flow = s_flow - 1_pI64
    if (flow(s_flow-1_pI64:s_flow-1_pI64) == ',') s_flow = s_flow - 1_pI64
    flow(s_flow:s_flow) = '}'
    s_flow = s_flow + 1_pI64
  else
    call line_toFlow(flow,s_flow,line)
  end if

end subroutine line_isFlow


!-------------------------------------------------------------------------------------------------
! @brief Transform a line of YAML of type <key>: <value> to flow style.
! @details Ensures that the <value> is consistent with the input required in the DAMASK YAML parser.
!-------------------------------------------------------------------------------------------------
recursive subroutine keyValue_toFlow(flow,s_flow,line)

  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer(pI64),    intent(inout) :: s_flow                                                         !< start position in flow
  character(len=*), intent(in)    :: line

  character(len=:), allocatable   :: line_asStandard                                                ! standard form of <key>: <value>
  integer(pI64) :: &
    d_flow, &
    col_pos, &
    offset_value


  col_pos = index(line,':',kind=pI64)
  if (line(col_pos+1_pI64:col_pos+1_pI64) /= ' ') call IO_error(704,ext_msg=line)
  if (isFlow(line(col_pos+1_pI64:))) then
    d_flow = len_trim(adjustl(line(:col_pos)),pI64)
    flow(s_flow:s_flow+d_flow+1_pI64) = trim(adjustl(line(:col_pos)))//' '
    s_flow = s_flow + d_flow + 1_pI64
    call line_isFlow(flow,s_flow,line(col_pos+1_pI64:))
  else
    offset_value = indentDepth(line(col_pos+2_pI64:))
    line_asStandard = line(:col_pos+1_pI64)//line(col_pos+2_pI64+offset_value:)
    call line_toFlow(flow,s_flow,line_asStandard)
  end if

end subroutine keyValue_toFlow


!-------------------------------------------------------------------------------------------------
! @brief Transform a line of YAML to flow style.
!-------------------------------------------------------------------------------------------------
subroutine line_toFlow(flow,s_flow,line)

  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer(pI64),    intent(inout) :: s_flow                                                         !< start position in flow
  character(len=*), intent(in)    :: line

  integer(pI64) :: d_flow


  d_flow = len_trim(adjustl(line),pI64)
  flow(s_flow:s_flow+d_flow) = trim(adjustl(line))
  s_flow = s_flow + d_flow

end subroutine line_toFlow


!-------------------------------------------------------------------------------------------------
! @brief Transform a block-style list to flow style.
! @details enters the function when encountered with the list indicator '- '
! reads each scalar list item and separates each other with a ','
! If list item is non scalar, it stores the offset for that list item block
! Call the 'decide' function if there is an increase in the indentation level or the list item is not a scalar
! decrease in indentation level indicates the end of an indentation block
!-------------------------------------------------------------------------------------------------
recursive subroutine lst(blck,flow,s_blck,s_flow,offset)

  character(len=*), intent(in)    :: blck                                                           !< YAML in mixed style
  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer(pI64),    intent(inout) :: s_blck, &                                                      !< start position in blck
                                     s_flow, &                                                      !< start position in flow
                                     offset                                                         !< stores leading '- ' in nested lists
  character(len=:), allocatable :: line,flow_line,inline
  integer(pI64) :: e_blck,indent


  indent = indentDepth(blck(s_blck:),offset)
  do while (s_blck <= len_trim(blck,pI64))
    e_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64
    line = clean(blck(s_blck:e_blck))
    if (trim(line) == '---' .or. trim(line) == '...') then
      exit
    elseif (len_trim(line,pI64) == 0_pI64) then
      s_blck = e_blck + 2_pI64                                                                      ! forward to next line
      cycle
    elseif (indentDepth(line,offset) > indent) then
      call decide(blck,flow,s_blck,s_flow,offset)
      offset = 0_pI64
      flow(s_flow:s_flow+1) = ', '
      s_flow = s_flow + 2_pI64
    elseif (indentDepth(line,offset) < indent .or. .not. isListItem(line)) then
      offset = 0_pI64
      exit                                                                                          ! job done (lower level)
    else
      if (trim(adjustl(line)) == '-') then                                                          ! list item in next line
        s_blck = e_blck + 2_pI64
        call skip_empty_lines(blck,s_blck)
        e_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64
        line = clean(blck(s_blck:e_blck))
        if (trim(line) == '---') call IO_error(707,ext_msg=line)
        if (indentDepth(line) < indent .or. indentDepth(line) == indent) &
          call IO_error(701,ext_msg=line)

        if (isScalar(line)) then
          call line_toFlow(flow,s_flow,line)
          s_blck = e_blck + 2_pI64
          offset = 0_pI64
        elseif (isFlow(line)) then
          if (isFlowList(line)) then
            call remove_line_break(blck,s_blck,']',flow_line)
          else
            call remove_line_break(blck,s_blck,'}',flow_line)
          end if
          call line_isFlow(flow,s_flow,flow_line)
          offset = 0_pI64
        end if
      else                                                                                          ! list item in the same line
        line = line(indentDepth(line)+3_pI64:)
        if (isScalar(line)) then
          call list_item_inline(blck,s_blck,inline,offset)
          offset = 0_pI64
          call line_toFlow(flow,s_flow,inline)
        elseif (isFlow(line)) then
          s_blck = s_blck + index(blck(s_blck:),'-',kind=pI64)
          if (isFlowList(line)) then
            call remove_line_break(blck,s_blck,']',flow_line)
          else
            call remove_line_break(blck,s_blck,'}',flow_line)
          end if
          call line_isFlow(flow,s_flow,flow_line)
          offset = 0_pI64
        else                                                                                        ! non scalar list item
          offset = offset + indentDepth(blck(s_blck:))+1_pI64                                       ! offset in spaces to be ignored
          s_blck = s_blck + index(blck(s_blck:e_blck),'-',kind=pI64)                                ! s_blck after '-' symbol
         end if
      end if
    end if

    if (isScalar(line) .or. isFlow(line)) then
      flow(s_flow:s_flow+1_pI64) = ', '
      s_flow = s_flow + 2_pI64
    end if

  end do

  s_flow = s_flow-1_pI64
  if (flow(s_flow-1_pI64:s_flow-1_pI64) == ',') s_flow = s_flow-1_pI64

end subroutine lst


!--------------------------------------------------------------------------------------------------
! @brief Transform a block-style dict to flow style.
! @details enters the function when encountered with the dictionary indicator ':'
! parses each line in the block and compares indentation of a line with the preceding line
! upon increase in indentation level -> 'decide' function decides if the line is a list or dict
! decrease in indentation indicates the end of an indentation block
!--------------------------------------------------------------------------------------------------
recursive subroutine dct(blck,flow,s_blck,s_flow,offset)

  character(len=*), intent(in)    :: blck                                                           !< YAML in mixed style
  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer(pI64),    intent(inout) :: s_blck, &                                                      !< start position in blck
                                     s_flow, &                                                      !< start position in flow
                                     offset

  character(len=:), allocatable :: line,flow_line
  integer(pI64) :: e_blck,indent,col_pos
  logical :: previous_isKey

  previous_isKey = .false.

  indent = indentDepth(blck(s_blck:),offset)

  do while (s_blck <= len_trim(blck,pI64))
    e_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64
    line = clean(blck(s_blck:e_blck))
    if (trim(line) == '---' .or. trim(line) == '...') then
      exit
    elseif (len_trim(line,pI64) == 0_pI64) then
      s_blck = e_blck + 2_pI64                                                                      ! forward to next line
      cycle
    elseif (indentDepth(line,offset) < indent) then
      if (isScalar(line) .or. isFlow(line) .and. previous_isKey) &
        call IO_error(701,ext_msg=line)
      offset = 0_pI64
      exit                                                                                          ! job done (lower level)
    elseif (indentDepth(line,offset) > indent .or. isListItem(line)) then
      offset = 0_pI64
      call decide(blck,flow,s_blck,s_flow,offset)
    else
      if (isScalar(line)) call IO_error(701,ext_msg=line)
      if (isFlow(line))   call IO_error(702,ext_msg=line)

      line = line(indentDepth(line)+1_pI64:)
      if (previous_isKey) then
        flow(s_flow-1_pI64:s_flow) = ', '
        s_flow = s_flow + 1_pI64
      end if

      if (isKeyValue(line)) then
        col_pos = index(line,':',kind=pI64)
        if (isFlow(line(col_pos+1_pI64:))) then
          if (isFlowList(line(col_pos+1_pI64:))) then
            call remove_line_break(blck,s_blck,']',flow_line)
          else
            call remove_line_break(blck,s_blck,'}',flow_line)
          end if
          call keyValue_toFlow(flow,s_flow,flow_line)
        else
          call keyValue_toFlow(flow,s_flow,line)
          s_blck = e_blck + 2_pI64
        end if
      else
        call line_toFlow(flow,s_flow,line)
        s_blck = e_blck + 2_pI64
      end if
    end if

    if (isScalar(line) .or. isKeyValue(line)) then
      flow(s_flow:s_flow) = ','
      s_flow = s_flow + 1_pI64
      previous_isKey = .false.
    else
      previous_isKey = .true.
    end if

    flow(s_flow:s_flow) = ' '
    s_flow = s_flow + 1_pI64
    offset = 0_pI64
  end do

  s_flow = s_flow - 1_pI64
  if (flow(s_flow-1_pI64:s_flow-1_pI64) == ',') s_flow = s_flow - 1_pI64

end subroutine dct


!--------------------------------------------------------------------------------------------------
! @brief Decide whether next block is list or dict.
!--------------------------------------------------------------------------------------------------
recursive subroutine decide(blck,flow,s_blck,s_flow,offset)

  character(len=*), intent(in)    :: blck                                                           !< YAML in mixed style
  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer(pI64),    intent(inout) :: s_blck, &                                                      !< start position in blck
                                     s_flow, &                                                      !< start position in flow
                                     offset
  integer(pI64) :: e_blck
  character(len=:), allocatable :: line,flow_line


  if (s_blck <= len(blck,kind=pI64)) then
    call skip_empty_lines(blck,s_blck)
    e_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64
    line = clean(blck(s_blck:e_blck))
    if (trim(line) == '---' .or. trim(line) == '...') then
      continue                                                                                      ! end parsing at this point but not stop the simulation
    elseif (len_trim(line,pI64) == 0_pI64) then
      s_blck = e_blck + 2_pI64
      call decide(blck,flow,s_blck,s_flow,offset)
    elseif (isListItem(line)) then
      flow(s_flow:s_flow) = '['
      s_flow = s_flow + 1_pI64
      call lst(blck,flow,s_blck,s_flow,offset)
      flow(s_flow:s_flow) = ']'
      s_flow = s_flow + 1_pI64
    elseif (isKey(line) .or. isKeyValue(line)) then
      flow(s_flow:s_flow) = '{'
      s_flow = s_flow + 1_pI64
      call dct(blck,flow,s_blck,s_flow,offset)
      flow(s_flow:s_flow) = '}'
      s_flow = s_flow + 1_pI64
    elseif (isFlow(line)) then
      if (isFlowList(line)) then
        call remove_line_break(blck,s_blck,']',flow_line)
      else
        call remove_line_break(blck,s_blck,'}',flow_line)
      end if
      call line_isFlow(flow,s_flow,line)
    else
      line = line(indentDepth(line)+1:)
      call line_toFlow(flow,s_flow,line)
      s_blck = e_blck + 2_pI64
    end if
  end if

end subroutine decide


!--------------------------------------------------------------------------------------------------
!> @brief Convert all block-style parts to flow style.
!> @details The input needs to end with a newline.
!--------------------------------------------------------------------------------------------------
function to_flow(blck)

  character(len=:), allocatable :: to_flow
  character(len=*), intent(in)  :: blck                                                             !< YAML mixed style

  character(len=:), allocatable :: line
  integer(pI64)                 :: s_blck, &                                                        !< start position in blck
                                   s_flow, &                                                        !< start position in flow
                                   offset, &                                                        !< counts leading '- ' in nested lists
                                   end_line

  allocate(character(len=len(blck,kind=pI64)*2)::to_flow)
  s_flow = 1_pI64
  s_blck = 1_pI64
  offset = 0_pI64

  if (len_trim(blck,pI64) /= 0_pI64) then
    call skip_empty_lines(blck,s_blck)
    call skip_file_header(blck,s_blck)
    line = clean(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL,kind=pI64) - 2_pI64))
    if (trim(line) == '---') s_blck = s_blck + index(blck(s_blck:),IO_EOL,kind=pI64)
    call decide(blck,to_flow,s_blck,s_flow,offset)
  end if
  line = clean(blck(s_blck:s_blck+index(blck(s_blck:),IO_EOL,kind=pI64)-2_pI64))
  if (trim(line)== '---') call IO_warning(709,ext_msg=line)
  to_flow = trim(to_flow(:s_flow-1_pI64))
  end_line = index(to_flow,IO_EOL,kind=pI64)
  if (end_line > 0_pI64) to_flow = to_flow(:end_line-1_pI64)

end function to_flow


!--------------------------------------------------------------------------------------------------
! @brief Remove comments (characters beyond '#') and trailing space.
!--------------------------------------------------------------------------------------------------
function clean(line)

  character(len=*), intent(in)  :: line
  character(len=:), allocatable :: clean

  integer(pI64) :: split
  character, parameter :: COMMENT_CHAR = '#'


  split = index(line,COMMENT_CHAR,kind=pI64)

  if (split == 0_pI64) then
    clean = trim(line)
  else
    clean = trim(line(:split-1_pI64))
  end if

end function clean


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some YAML functions.
!--------------------------------------------------------------------------------------------------
subroutine YAML_selfTest()

  if (indentDepth(' a') /= 1)     error stop 'indentDepth'
  if (indentDepth('a')  /= 0)     error stop 'indentDepth'
  if (indentDepth('x ') /= 0)     error stop 'indentDepth'

  if (.not. quotedStr("'a'"))     error stop 'quotedStr'

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

  if (       isScalar('a:  '))     error stop 'isScalar'
  if (       isScalar('a: b'))     error stop 'isScalar'
  if (       isScalar('{a:b}'))    error stop 'isScalar'
  if (       isScalar('- a:'))     error stop 'isScalar'
  if (.not.  isScalar('   a'))     error stop 'isScalar'

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
      '%YAML 1.1'//IO_EOL//&
      '---'//IO_EOL//&
      'a:     ["b",'//IO_EOL//&
      'c: '//IO_EOL//&
      '"d",                               "e"]'//IO_EOL

    character(len=*), parameter :: flow = &
      '{a: ["b", {c: "d"}, "e"]}'

    if ( .not. to_flow(flow_multi)        == flow) error stop 'to_flow'
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

    if ( .not. to_flow(flow_multi)        == flow) error stop 'to_flow'
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
      " - c:d"//IO_EOL//&
      "  e.f,"//IO_EOL//&
      " bb:"//IO_EOL//&
      " "//IO_EOL//&
      "  - "//IO_EOL//&
      "   {param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}"//IO_EOL//&
      "..."//IO_EOL
    character(len=*), parameter :: mixed_flow = &
      '{aa: [{param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}, "c:d e.f,"], bb: [{param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}]}'

    if (.not. to_flow(block_flow) == mixed_flow)    error stop 'to_flow'
  end block basic_mixed

  parse: block

    type(tDict), pointer :: dict
    type(tList), pointer :: list
    character(len=*), parameter :: &
      lst = '[1, 2, 3, 4]', &
      dct = '{a: 1, b: 2}'

    list => YAML_str_asList(lst//IO_EOL)
    if (list%asFormattedStr() /= lst) error stop 'str_asList'
    dict => YAML_str_asDict(dct//IO_EOL)
    if (dict%asFormattedStr() /= dct) error stop 'str_asDict'

  end block parse

  comment: block
    character(len=:),      allocatable :: str,out

    str='#';out=clean(str)
    if (out /= ''   .or. len(out) /= 0)  error stop 'clean/1'
    str=' #';out=clean(str)
    if (out /= ''   .or. len(out) /= 0)  error stop 'clean/2'
    str=' # ';out=clean(str)
    if (out /= ''   .or. len(out) /= 0)  error stop 'clean/3'
    str=' # a';out=clean(str)
    if (out /= ''   .or. len(out) /= 0)  error stop 'clean/4'
    str=' a#';out=clean(str)
    if (out /= ' a' .or. len(out) /= 2)  error stop 'clean/5'
    str=' ab #';out=clean(str)
    if (out /= ' ab'.or. len(out) /= 3)  error stop 'clean/6'
  end block comment

end subroutine YAML_selfTest
#endif

end module YAML
