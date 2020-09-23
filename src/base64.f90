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

  print'(/,a)', ' <<<+-  base64 init  -+>>>'; flush(IO_STDOUT)

  call selfTest

end subroutine base64_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate number of Base64 characters required for storage of N bytes.
!--------------------------------------------------------------------------------------------------
pure function base64_nChar(nByte)

  integer(pI64), intent(in) :: nByte
  integer(pI64)             :: base64_nChar

  base64_nChar = 4_pI64 * (nByte/3_pI64 + merge(1_pI64,0_pI64,mod(nByte,3_pI64) /= 0_pI64))

end function base64_nChar


!--------------------------------------------------------------------------------------------------
!> @brief Calculate number of bytes required for storage of N Base64 characters.
!--------------------------------------------------------------------------------------------------
pure function base64_nByte(nBase64)

  integer(pI64), intent(in) :: nBase64
  integer(pI64)             :: base64_nByte

  base64_nByte = 3_pI64 * (nBase64/4_pI64)

end function base64_nByte


!--------------------------------------------------------------------------------------------------
!> @brief Decode Base64 ASCII string into byte-wise binary representation.
!--------------------------------------------------------------------------------------------------
function base64_to_bytes(base64_str,s,e) result(bytes)

  character(len=*),  intent(in) :: base64_str                                                       !< Base64 string representation
  integer(pI64), intent(in), optional :: &
    s, &                                                                                            !< start (in bytes)
    e                                                                                               !< end (in bytes)

  integer(pI64) :: s_bytes, e_bytes, s_str, e_str
  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes

  if(.not. validBase64(base64_str)) call IO_error(114,ext_msg='invalid character')

  if(present(s)) then
    if(s<1_pI64) call IO_error(114, ext_msg='s out of range')
    s_str = ((s-1_pI64)/3_pI64)*4_pI64 + 1_pI64
    s_bytes = mod(s-1_pI64,3_pI64) + 1_pI64
  else
    s_str = 1_pI64
    s_bytes = 1_pI64
  endif

  if(present(e)) then
    if(e>base64_nByte(len(base64_str,kind=pI64))) call IO_error(114, ext_msg='e out of range')
    e_str = ((e-1_pI64)/3_pI64)*4_pI64 + 4_pI64
    e_bytes = e - base64_nByte(s_str)
  else
    e_str = len(base64_str,kind=pI64)
    e_bytes = base64_nByte(len(base64_str,kind=pI64)) - base64_nByte(s_str)
    if(base64_str(e_str-0_pI64:e_str-0_pI64) == '=') e_bytes = e_bytes - 1_pI64
    if(base64_str(e_str-1_pI64:e_str-1_pI64) == '=') e_bytes = e_bytes - 1_pI64
  endif

  bytes = decodeBase64(base64_str(s_str:e_str))
  bytes = bytes(s_bytes:e_bytes)

end function base64_to_bytes


!--------------------------------------------------------------------------------------------------
!> @brief Convert a Base64 ASCII string into its byte-wise binary representation.
!--------------------------------------------------------------------------------------------------
pure function decodeBase64(base64_str) result(bytes)

  character(len=*), intent(in) :: base64_str                                                        !< Base64 string representation

  integer(C_SIGNED_CHAR), dimension(base64_nByte(len(base64_str,pI64))) :: bytes

  integer(C_SIGNED_CHAR), dimension(0:3) :: charPos
  integer(pI64) :: c, b, p

  c = 1_pI64
  b = 1_pI64

  do while(c < len(base64_str,kind=pI64))
    do p=0_pI64,3_pI64
      if(c+p<=len(base64_str,kind=pI64)) then
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
    b = b+3_pI64
    c = c+4_pI64
  enddo

end function decodeBase64


!--------------------------------------------------------------------------------------------------
!> @brief Test for valid Base64 encoded string.
!> @details Input string must be properly padded.
!--------------------------------------------------------------------------------------------------
pure logical function validBase64(base64_str)

  character(len=*), intent(in) :: base64_str                                                        !< Base64 string representation

  integer(pI64) :: l

  l = len(base64_str,pI64)
  validBase64 = .true.

  if(mod(l,4_pI64)/=0_pI64 .or. l < 4_pInt)                                  validBase64 = .false.
  if(verify(base64_str(:l-2_pI64),base64_encoding,     kind=pI64) /= 0_pI64) validBase64 = .false.
  if(verify(base64_str(l-1_pI64:),base64_encoding//'=',kind=pI64) /= 0_pI64) validBase64 = .false.

end function validBase64


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of base64 functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes
  character(len=*), parameter :: zero_to_three = 'AAECAw=='

  ! https://en.wikipedia.org/wiki/Base64#Output_padding
  if(base64_nChar(20_pI64) /= 28_pI64) error stop 'base64_nChar/20/28'
  if(base64_nChar(19_pI64) /= 28_pI64) error stop 'base64_nChar/19/28'
  if(base64_nChar(18_pI64) /= 24_pI64) error stop 'base64_nChar/18/24'
  if(base64_nChar(17_pI64) /= 24_pI64) error stop 'base64_nChar/17/24'
  if(base64_nChar(16_pI64) /= 24_pI64) error stop 'base64_nChar/16/24'

  if(base64_nByte(4_pI64)  /= 3_pI64)  error stop 'base64_nByte/4/3'
  if(base64_nByte(8_pI64)  /= 6_pI64)  error stop 'base64_nByte/8/6'

  bytes = base64_to_bytes(zero_to_three)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) error stop 'base64_to_bytes//'

  bytes = base64_to_bytes(zero_to_three,e=1_pI64)
  if(any(bytes /= int([0],C_SIGNED_CHAR))       .or. size(bytes) /= 1) error stop 'base64_to_bytes//1'
  bytes = base64_to_bytes(zero_to_three,e=2_pI64)
  if(any(bytes /= int([0,1],C_SIGNED_CHAR))     .or. size(bytes) /= 2) error stop 'base64_to_bytes//2'
  bytes = base64_to_bytes(zero_to_three,e=3_pI64)
  if(any(bytes /= int([0,1,2],C_SIGNED_CHAR))   .or. size(bytes) /= 3) error stop 'base64_to_bytes//3'
  bytes = base64_to_bytes(zero_to_three,e=4_pI64)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) error stop 'base64_to_bytes//4'

  bytes = base64_to_bytes(zero_to_three,s=1_pI64)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) error stop 'base64_to_bytes/1/'
  bytes = base64_to_bytes(zero_to_three,s=2_pI64)
  if(any(bytes /= int([1,2,3],C_SIGNED_CHAR))   .or. size(bytes) /= 3) error stop 'base64_to_bytes/2/'
  bytes = base64_to_bytes(zero_to_three,s=3_pI64)
  if(any(bytes /= int([2,3],C_SIGNED_CHAR))     .or. size(bytes) /= 2) error stop 'base64_to_bytes/3/'
  bytes = base64_to_bytes(zero_to_three,s=4_pI64)
  if(any(bytes /= int([3],C_SIGNED_CHAR))       .or. size(bytes) /= 1) error stop 'base64_to_bytes/4/'

  bytes = base64_to_bytes(zero_to_three,s=1_pI64,e=1_pI64)
  if(any(bytes /= int([0],C_SIGNED_CHAR))       .or. size(bytes) /= 1) error stop 'base64_to_bytes/1/1'
  bytes = base64_to_bytes(zero_to_three,s=2_pI64,e=2_pI64)
  if(any(bytes /= int([1],C_SIGNED_CHAR))       .or. size(bytes) /= 1) error stop 'base64_to_bytes/2/2'
  bytes = base64_to_bytes(zero_to_three,s=3_pI64,e=3_pI64)
  if(any(bytes /= int([2],C_SIGNED_CHAR))       .or. size(bytes) /= 1) error stop 'base64_to_bytes/3/3'
  bytes = base64_to_bytes(zero_to_three,s=4_pI64,e=4_pI64)
  if(any(bytes /= int([3],C_SIGNED_CHAR))       .or. size(bytes) /= 1) error stop 'base64_to_bytes/4/4'

  bytes = base64_to_bytes(zero_to_three,s=1_pI64,e=2_pI64)
  if(any(bytes /= int([0,1],C_SIGNED_CHAR))     .or. size(bytes) /= 2) error stop 'base64_to_bytes/1/2'
  bytes = base64_to_bytes(zero_to_three,s=2_pI64,e=3_pI64)
  if(any(bytes /= int([1,2],C_SIGNED_CHAR))     .or. size(bytes) /= 2) error stop 'base64_to_bytes/2/3'
  bytes = base64_to_bytes(zero_to_three,s=3_pI64,e=4_pI64)
  if(any(bytes /= int([2,3],C_SIGNED_CHAR))     .or. size(bytes) /= 2) error stop 'base64_to_bytes/3/4'

  bytes = base64_to_bytes(zero_to_three,s=1_pI64,e=3_pI64)
  if(any(bytes /= int([0,1,2],C_SIGNED_CHAR))   .or. size(bytes) /= 3) error stop 'base64_to_bytes/1/3'
  bytes = base64_to_bytes(zero_to_three,s=2_pI64,e=4_pI64)
  if(any(bytes /= int([1,2,3],C_SIGNED_CHAR))   .or. size(bytes) /= 3) error stop 'base64_to_bytes/2/4'

  bytes = base64_to_bytes(zero_to_three,s=1_pI64,e=4_pI64)
  if(any(bytes /= int([0,1,2,3],C_SIGNED_CHAR)) .or. size(bytes) /= 4) error stop 'base64_to_bytes/1/4'

end subroutine selfTest

end module base64
