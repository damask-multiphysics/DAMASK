!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Decode Base64 strings.
!> @details See https://en.wikipedia.org/wiki/Base64.
!--------------------------------------------------------------------------------------------------
module base64
  use prec
  use IO

  implicit none
  private

  character(len=*), parameter :: &
    base64_encoding='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'

  public :: &
    base64_init, &
    base64_to_bytes, &
    base64_nChar, &
    base64_nByte

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine base64_init

  write(6,'(/,a)') ' <<<+-  base64 init  -+>>>'; flush(6)

  call selfTest

end subroutine base64_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate number of Base64 characters required for storage of N bytes.
!--------------------------------------------------------------------------------------------------
pure function base64_nChar(nByte)

  integer(pLongInt), intent(in) :: nByte
  integer(pLongInt)             :: base64_nChar

  base64_nChar = 4_pLongInt * (nByte/3_pLongInt + merge(1_pLongInt,0_pLongInt,mod(nByte,3_pLongInt) /= 0_pLongInt))

end function base64_nChar


!--------------------------------------------------------------------------------------------------
!> @brief Calculate number of bytes required for storage of N Base64 characters.
!--------------------------------------------------------------------------------------------------
pure function base64_nByte(nBase64)

  integer(pLongInt), intent(in) :: nBase64
  integer(pLongInt)             :: base64_nByte

  base64_nByte = 3_pLongInt * (nBase64/4_pLongInt)

end function base64_nByte


!--------------------------------------------------------------------------------------------------
!> @brief Decode Base64 ASCII string into byte-wise binary representation.
!--------------------------------------------------------------------------------------------------
function base64_to_bytes(base64_str,s,e) result(bytes)

  character(len=*),  intent(in) :: base64_str                                                       !< Base64 string representation
  integer(pLongInt), intent(in), optional :: &
    s, &                                                                                            !< start (in bytes)
    e                                                                                               !< end (in bytes)

  integer(pLongInt) :: s_bytes, e_bytes, s_str, e_str
  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes

  if(.not. valid_base64(base64_str)) call IO_error(114,ext_msg='invalid character')

  if(present(s)) then
    if(s<1_pLongInt) call IO_error(114, ext_msg='s out of range')
    s_str = ((s-1_pLongInt)/3_pLongInt)*4_pLongInt + 1_pLongInt
    s_bytes = mod(s-1_pLongInt,3_pLongInt) + 1_pLongInt
  else
    s_str = 1_pLongInt
    s_bytes = 1_pLongInt
  endif

  if(present(e)) then
    if(e>base64_nByte(len(base64_str,kind=pLongInt))) call IO_error(114, ext_msg='e out of range')
    e_str = ((e-1_pLongInt)/3_pLongInt)*4_pLongInt + 4_pLongInt
    e_bytes = e - base64_nByte(s_str)
  else
    e_str = len(base64_str,kind=pLongInt)
    e_bytes = base64_nByte(len(base64_str,kind=pLongInt)) - base64_nByte(s_str)
    if(base64_str(e_str-0_pLongInt:e_str-0_pLongInt) == '=') e_bytes = e_bytes - 1_pLongInt
    if(base64_str(e_str-1_pLongInt:e_str-1_pLongInt) == '=') e_bytes = e_bytes - 1_pLongInt
  endif

  bytes = decode_base64(base64_str(s_str:e_str))
  bytes = bytes(s_bytes:e_bytes)

end function base64_to_bytes


!--------------------------------------------------------------------------------------------------
!> @brief Convert a Base64 ASCII string into its byte-wise binary representation.
!--------------------------------------------------------------------------------------------------
pure function decode_base64(base64_str) result(bytes)

  character(len=*), intent(in) :: base64_str                                                        !< Base64 string representation

  integer(C_SIGNED_CHAR), dimension(base64_nByte(len(base64_str,pLongInt))) :: bytes

  integer(C_SIGNED_CHAR), dimension(0:3) :: charPos
  integer(pLongInt) :: c, b, p

  c = 1_pLongInt
  b = 1_pLongInt

  do while(c < len(base64_str,kind=pLongInt))
    do p=0_pLongInt,3_pLongInt
      if(c+p<=len(base64_str,kind=pLongInt)) then
        charPos(p) = int(index(base64_encoding,base64_str(c+p:c+p))-1,C_SIGNED_CHAR)
      else
        charPos(p) = 0_C_SIGNED_CHAR
      endif
    enddo

    call mvbits(charPos(0),0,6,bytes(b+0),2)
    call mvbits(charPos(1),4,2,bytes(b+0),0)
    call mvbits(charPos(1),0,4,bytes(b+1),4)
    call mvbits(charPos(2),2,4,bytes(b+1),0)
    call mvbits(charPos(2),0,2,bytes(b+2),6)
    call mvbits(charPos(3),0,6,bytes(b+2),0)
    b = b+3_pLongInt
    c = c+4_pLongInt
  enddo

end function decode_base64


!--------------------------------------------------------------------------------------------------
!> @brief Test for valid Base64 encoded string.
!> @details Input string must be properly padded.
!--------------------------------------------------------------------------------------------------
pure logical function valid_base64(base64_str)

  character(len=*), intent(in) :: base64_str                                                        !< Base64 string representation

  integer(pLongInt) :: l

  l = len(base64_str,pLongInt)
  valid_base64 = .true.

  if(mod(l,4_pLongInt)/=0_pLongInt .or. l < 4_pInt)                                      valid_base64 = .false.
  if(verify(base64_str(:l-2_pLongInt),base64_encoding,     kind=pLongInt) /= 0_pLongInt) valid_base64 = .false.
  if(verify(base64_str(l-1_pLongInt:),base64_encoding//'=',kind=pLongInt) /= 0_pLongInt) valid_base64 = .false.

end function valid_base64


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of base64 functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes
  character(len=*), parameter :: zero_to_three = 'AAECAw=='

  ! https://en.wikipedia.org/wiki/Base64#Output_padding
  if(base64_nChar(20_pLongInt) /= 28_pLongInt) call IO_error(0,ext_msg='base64_nChar/20/28')
  if(base64_nChar(19_pLongInt) /= 28_pLongInt) call IO_error(0,ext_msg='base64_nChar/19/28')
  if(base64_nChar(18_pLongInt) /= 24_pLongInt) call IO_error(0,ext_msg='base64_nChar/18/24')
  if(base64_nChar(17_pLongInt) /= 24_pLongInt) call IO_error(0,ext_msg='base64_nChar/17/24')
  if(base64_nChar(16_pLongInt) /= 24_pLongInt) call IO_error(0,ext_msg='base64_nChar/16/24')

  if(base64_nByte(4_pLongInt)  /= 3_pLongInt)  call IO_error(0,ext_msg='base64_nByte/4/3')
  if(base64_nByte(8_pLongInt)  /= 6_pLongInt)  call IO_error(0,ext_msg='base64_nByte/8/6')

  bytes = base64_to_bytes(zero_to_three)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) call IO_error(0,ext_msg='base64_to_bytes//')

  bytes = base64_to_bytes(zero_to_three,e=1_pLongInt)
  if(any(bytes /= int([0],C_SIGNED_CHAR))       .or. size(bytes) /= 1) call IO_error(0,ext_msg='base64_to_bytes//1')
  bytes = base64_to_bytes(zero_to_three,e=2_pLongInt)
  if(any(bytes /= int([0,1],C_SIGNED_CHAR))     .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes//2')
  bytes = base64_to_bytes(zero_to_three,e=3_pLongInt)
  if(any(bytes /= int([0,1,2],C_SIGNED_CHAR))   .or. size(bytes) /= 3) call IO_error(0,ext_msg='base64_to_bytes//3')
  bytes = base64_to_bytes(zero_to_three,e=4_pLongInt)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) call IO_error(0,ext_msg='base64_to_bytes//4')

  bytes = base64_to_bytes(zero_to_three,s=1_pLongInt)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) call IO_error(0,ext_msg='base64_to_bytes/1/')
  bytes = base64_to_bytes(zero_to_three,s=2_pLongInt)
  if(any(bytes /= int([1,2,3],C_SIGNED_CHAR))   .or. size(bytes) /= 3) call IO_error(0,ext_msg='base64_to_bytes/2/')
  bytes = base64_to_bytes(zero_to_three,s=3_pLongInt)
  if(any(bytes /= int([2,3],C_SIGNED_CHAR))     .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes/3/')
  bytes = base64_to_bytes(zero_to_three,s=4_pLongInt)
  if(any(bytes /= int([3],C_SIGNED_CHAR))       .or. size(bytes) /= 1) call IO_error(0,ext_msg='base64_to_bytes/4/')

  bytes = base64_to_bytes(zero_to_three,s=1_pLongInt,e=1_pLongInt)
  if(any(bytes /= int([0],C_SIGNED_CHAR))       .or. size(bytes) /= 1) call IO_error(0,ext_msg='base64_to_bytes/1/1')
  bytes = base64_to_bytes(zero_to_three,s=2_pLongInt,e=2_pLongInt)
  if(any(bytes /= int([1],C_SIGNED_CHAR))       .or. size(bytes) /= 1) call IO_error(0,ext_msg='base64_to_bytes/2/2')
  bytes = base64_to_bytes(zero_to_three,s=3_pLongInt,e=3_pLongInt)
  if(any(bytes /= int([2],C_SIGNED_CHAR))       .or. size(bytes) /= 1) call IO_error(0,ext_msg='base64_to_bytes/3/3')
  bytes = base64_to_bytes(zero_to_three,s=4_pLongInt,e=4_pLongInt)
  if(any(bytes /= int([3],C_SIGNED_CHAR))       .or. size(bytes) /= 1) call IO_error(0,ext_msg='base64_to_bytes/4/4')

  bytes = base64_to_bytes(zero_to_three,s=1_pLongInt,e=2_pLongInt)
  if(any(bytes /= int([0,1],C_SIGNED_CHAR))     .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes/1/2')
  bytes = base64_to_bytes(zero_to_three,s=2_pLongInt,e=3_pLongInt)
  if(any(bytes /= int([1,2],C_SIGNED_CHAR))     .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes/2/3')
  bytes = base64_to_bytes(zero_to_three,s=3_pLongInt,e=4_pLongInt)
  if(any(bytes /= int([2,3],C_SIGNED_CHAR))     .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes/3/4')

  bytes = base64_to_bytes(zero_to_three,s=1_pLongInt,e=3_pLongInt)
  if(any(bytes /= int([0,1,2],C_SIGNED_CHAR))   .or. size(bytes) /= 3) call IO_error(0,ext_msg='base64_to_bytes/1/3')
  bytes = base64_to_bytes(zero_to_three,s=2_pLongInt,e=4_pLongInt)
  if(any(bytes /= int([1,2,3],C_SIGNED_CHAR))   .or. size(bytes) /= 3) call IO_error(0,ext_msg='base64_to_bytes/2/4')

  bytes = base64_to_bytes(zero_to_three,s=1_pLongInt,e=4_pLongInt)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) call IO_error(0,ext_msg='base64_to_bytes/1/4')

end subroutine selfTest

end module base64
