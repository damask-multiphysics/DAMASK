!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Decode Base64 strings
!--------------------------------------------------------------------------------------------------
module base64
  use prec
  use IO

  implicit none
  private

  public :: &
    base64_init, &
    base64_to_bytes, &
    base64_nBase64, &
    base64_nByte

contains


!--------------------------------------------------------------------------------------------------
!> @brief do self test
!--------------------------------------------------------------------------------------------------
subroutine base64_init

  write(6,'(/,a)') ' <<<+-  base64 init  -+>>>'; flush(6)

  call selfTest

end subroutine base64_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate number of Base64 characters required for storage of N bytes
!--------------------------------------------------------------------------------------------------
pure function base64_nBase64(nByte)

  integer(pLongInt), intent(in) :: nByte
  integer(pLongInt)             :: base64_nBase64

  base64_nBase64 = 4_pLongInt * (nByte/3_pLongInt + merge(1_pLongInt,0_pLongInt,mod(nByte,3_pLongInt) /= 0_pLongInt))

end function base64_nBase64


!--------------------------------------------------------------------------------------------------
!> @brief calculate number of bytes required for storage of N Base64 characters
!--------------------------------------------------------------------------------------------------
pure function base64_nByte(nBase64)

  integer(pLongInt), intent(in) :: nBase64
  integer(pLongInt)             :: base64_nByte

  base64_nByte = 3_pLongInt * (nBase64/4_pLongInt)

end function base64_nByte


!--------------------------------------------------------------------------------------------------
!> @brief return byte-wise binary representation of a (part of a) Base64 ASCII string
!--------------------------------------------------------------------------------------------------
pure function base64_to_bytes(base64_str,s,e) result(bytes)

  character(len=*),  intent(in) :: base64_str                                                       !< Base64 string representation
  integer(pLongInt), intent(in), optional :: &
    s, &                                                                                            !< start (in bytes)
    e                                                                                               !< end (in bytes)

  integer(pLongInt) :: s_, e_
  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes, bytes_raw

  allocate(bytes_raw(base64_nByte(len_trim(base64_str,pLongInt))))

  if(present(s)) then
    s_ = s
  else
    s_ = 1
  endif

  if(present(e)) then
    e_ = e
  else
    e_ = size(bytes_raw,kind=pLongInt)
  endif

  ! ToDo: inefficient for s_>>1 or e_<<size(bytes_raw)
  call decode_base64(bytes_raw,base64_str)
  bytes = bytes_raw(s_:e_)

end function base64_to_bytes


!--------------------------------------------------------------------------------------------------
!> @brief convert a Base64 ASCII string into its byte-wise binary representation
!> @details https://en.wikipedia.org/wiki/Base64
!--------------------------------------------------------------------------------------------------
pure subroutine decode_base64(bytes,base64_str)

  integer(C_SIGNED_CHAR), intent(out),dimension(:) :: bytes                                         !< byte-wise representation
  character(len=*),       intent(in)               :: base64_str                                    !< Base64 string representation

  integer(C_SIGNED_CHAR), dimension(0:3) :: charPos
  integer(pLongInt) :: c, b, bytesLen,base64Len, p
  character(len=*), parameter :: encoding='ABCDEFGHIJKLMNOPQRSTUVWXYZ&
                                          &abcdefghijklmnopqrstuvwxyz0123456789+/'

  bytes = 0_C_SIGNED_CHAR

  bytesLen  = size(bytes,kind=pLongInt)
  base64Len = len_trim(base64_str,kind=pLongInt)

  c = 1_pLongInt
  b = 1_pLongInt

  do while(.True.)
    do p=0_pLongInt,3_pLongInt
      if(c+p<=len(base64_str,kind=pLongInt)) then
        charPos(p) = int(merge(index(encoding,base64_str(c+p:c+p))-1, 0, base64_str(c+p:c+p) /= '='),C_SIGNED_CHAR)
      else
        charPos(p) = 0_C_SIGNED_CHAR
      endif
    enddo

    if (b+0<=bytesLen) then
       call mvbits(charPos(0),0,6,bytes(b+0),2)
       call mvbits(charPos(1),4,2,bytes(b+0),0)
    endif
    if (b+1<=bytesLen) then
       call mvbits(charPos(1),0,4,bytes(b+1),4)
       call mvbits(charPos(2),2,4,bytes(b+1),0)
    endif
    if (b+2<=bytesLen) then
       call mvbits(charPos(2),0,2,bytes(b+2),6)
       call mvbits(charPos(3),0,6,bytes(b+2),0)
    endif
    b = b+3_pLongInt
    c = c+4_pLongInt
    if(b>bytesLen+3_pLongInt .or. c > base64Len+4_pLongInt) exit
  enddo

end subroutine decode_base64


!--------------------------------------------------------------------------------------------------
!> @brief check correctness of base64 functions
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes
  character(len=*), parameter :: zero_to_three = 'AAECAw=='

  allocate(bytes(7), source = -1_C_SIGNED_CHAR)
  call decode_base64(bytes,zero_to_three)
  if(any(bytes /= int([0,1,2,3,0,0,0],C_SIGNED_CHAR)) .or. size(bytes) /= 7) call IO_error(0,ext_msg='base64_decode')

  bytes = base64_to_bytes(zero_to_three)
  if(any(bytes /= int([0,1,2,3,0,0],C_SIGNED_CHAR))   .or. size(bytes) /= 6) call IO_error(0,ext_msg='base64_to_bytes')

  bytes = base64_to_bytes(zero_to_three,s=2_pLongInt)
  if(any(bytes /= int([1,2,3,0,0],C_SIGNED_CHAR))     .or. size(bytes) /= 5) call IO_error(0,ext_msg='base64_to_bytes/s')

  bytes = base64_to_bytes(zero_to_three,s=2_pLongInt,e=3_pLongInt)
  if(any(bytes /= int([1,2],C_SIGNED_CHAR))           .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes/s/e')

  bytes = base64_to_bytes(zero_to_three,e=2_pLongInt)
  if(any(bytes /= int([0,1],C_SIGNED_CHAR))           .or. size(bytes) /= 2) call IO_error(0,ext_msg='base64_to_bytes/e')

  ! https://en.wikipedia.org/wiki/Base64#Output_padding
  if(base64_nBase64(20_pLongInt) /= 28_pLongInt) call IO_error(0,ext_msg='base64_nBase64/20/28')
  if(base64_nBase64(19_pLongInt) /= 28_pLongInt) call IO_error(0,ext_msg='base64_nBase64/19/28')
  if(base64_nBase64(18_pLongInt) /= 24_pLongInt) call IO_error(0,ext_msg='base64_nBase64/18/24')
  if(base64_nBase64(17_pLongInt) /= 24_pLongInt) call IO_error(0,ext_msg='base64_nBase64/17/24')
  if(base64_nBase64(16_pLongInt) /= 24_pLongInt) call IO_error(0,ext_msg='base64_nBase64/16/24')

  if(base64_nByte(4_pLongInt)    /= 3_pLongInt)  call IO_error(0,ext_msg='base64_nByte/4/3')
  if(base64_nByte(8_pLongInt)    /= 6_pLongInt)  call IO_error(0,ext_msg='base64_nByte/8/6')

end subroutine selfTest

end module base64
