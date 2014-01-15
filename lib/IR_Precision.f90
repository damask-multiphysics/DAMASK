!> @addtogroup GlobalVarPar Global Variables and Parameters
!> List of global variables and parameters.
!> @addtogroup Interface Interfaces
!> List of explicitly defined interface.
!> @addtogroup Library Modules Libraries
!> List of modules containing libraries of procedures.
!> @addtogroup PublicProcedure Public Procedures
!> List of public procedures.
!> @addtogroup PrivateProcedure Private Procedures
!> List of private procedures.

!> @ingroup Library
!> @{
!> @defgroup IR_PrecisionLibrary IR_Precision
!> Portable kind-parameters module
!> @}

!> @ingroup Interface
!> @{
!> @defgroup IR_PrecisionInterface IR_Precision
!> Portable kind-parameters module
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup IR_PrecisionGlobalVarPar IR_Precision
!> Portable kind-parameters module
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup IR_PrecisionPublicProcedure IR_Precision
!> Portable kind-parameters module
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup IR_PrecisionPrivateProcedure IR_Precision
!> Portable kind-parameters module
!> @}

!> @brief     Module IR_Precision makes available some portable kind-parameters and some useful procedures to deal with them.
!> @details   It also provides variables that contain the minimum and maximum representable values, smallest real values and
!>            smallest representable differences by the running calculator.
!>
!>            Finally the module provides procedures to convert a string to number and vice versa, a function to check the endianism
!>            of the running calculator and a procedure to print all the aboves values.
!> @note The \em quadruple precision R16P could be activated defining \b r16p pre-processor flag (e.g. -Dr16p). Furthermore if
!> compiling with Portland Group Compiler define the pre-processor flag \b pgf95 to avoid error in computing \em Zero variables:
!> pgf compiler doesn't accept \b nearest built-in function in variables initialization.
!> @author    Stefano Zaghi
!> @version   1.0
!> @date      2012-04-24
!> @copyright GNU Public License version 3.
!> @todo \b g95_test: Test g95 compiler
!> @ingroup IR_PrecisionLibrary
module IR_Precision
!-----------------------------------------------------------------------------------------------------------------------------------
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: endianL,endianB,endian
public:: R16P, FR16P, DR16P, MinR16P, MaxR16P, BIR16P, BYR16P, smallR16P, ZeroR16
public:: R8P,  FR8P,  DR8P,  MinR8P,  MaxR8P,  BIR8P,  BYR8P,  smallR8P,  ZeroR8
public:: R4P,  FR4P,  DR4P,  MinR4P,  MaxR4P,  BIR4P,  BYR4P,  smallR4P,  ZeroR4
public:: R_P,  FR_P,  DR_P,  MinR_P,  MaxR_P,  BIR_P,  BYR_P,  smallR_P,  Zero
public:: I8P,  FI8P,  DI8P,  MinI8P,  MaxI8P,  BII8P,  BYI8P
public:: I4P,  FI4P,  DI4P,  MinI4P,  MaxI4P,  BII4P,  BYI4P
public:: I2P,  FI2P,  DI2P,  MinI2P,  MaxI2P,  BII2P,  BYI2P
public:: I1P,  FI1P,  DI1P,  MinI1P,  MaxI1P,  BII1P,  BYI1P
public:: I_P,  FI_P,  DI_P,  MinI_P,  MaxI_P,  BII_P,  BYI_P
public:: NRknd, RPl, FRl
public:: NIknd, RIl, FIl
public:: check_endian
public:: bit_size,byte_size
public:: str, strz, cton, bstr, bcton
public:: ir_initialized,IR_Init
public:: IR_Print
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup IR_PrecisionGlobalVarPar
!> @{
logical:: ir_initialized = .false. !< Flag for chcecking the initialization of some variables that must be initialized by IR_Init.
! Bit ordering of the running architecture:
integer, parameter:: endianL = 1       !< Little endian parameter.
integer, parameter:: endianB = 0       !< Big endian parameter.
integer::            endian  = endianL !< Bit ordering: Little endian (endianL), or Big endian (endianB).

! The following are the portable kind parameters available.
! Real precision definitions:
#ifdef r16p
integer, parameter:: R16P = selected_real_kind(33,4931) !< 33  digits, range \f$[10^{-4931}, 10^{+4931} - 1]\f$; 128 bits.
#else
integer, parameter:: R16P = selected_real_kind(15,307)  !< Defined as R8P; 64 bits.
#endif
integer, parameter:: R8P  = selected_real_kind(15,307)  !< 15  digits, range \f$[10^{-307} , 10^{+307}  - 1]\f$; 64 bits.
integer, parameter:: R4P  = selected_real_kind(6,37)    !< 6   digits, range \f$[10^{-37}  , 10^{+37}   - 1]\f$; 32 bits.
integer, parameter:: R_P  = R8P                         !< Default real precision.
! Integer precision definitions:
integer, parameter:: I8P  = selected_int_kind(18) !< Range \f$[-2^{63},+2^{63} - 1]\f$, 19 digits plus sign; 64 bits.
integer, parameter:: I4P  = selected_int_kind(9)  !< Range \f$[-2^{31},+2^{31} - 1]\f$, 10 digits plus sign; 32 bits.
integer, parameter:: I2P  = selected_int_kind(4)  !< Range \f$[-2^{15},+2^{15} - 1]\f$, 5  digits plus sign; 16 bits.
integer, parameter:: I1P  = selected_int_kind(2)  !< Range \f$[-2^{7} ,+2^{7}  - 1]\f$, 3  digits plus sign; 8  bits.
integer, parameter:: I_P  = I4P                   !< Default integer precision.

! Format parameters useful for writing in a well-ascii-format numeric variables.
! Real output formats:
character(10), parameter:: FR16P = '(E42.33E4)' !< Output format for kind=R16P variable.
character(10), parameter:: FR8P  = '(E23.15E3)' !< Output format for kind=R8P variable.
character(9),  parameter:: FR4P  = '(E13.6E2)'  !< Output format for kind=R4P variable.
character(10), parameter:: FR_P  = FR8P         !< Output format for kind=R_P variable.
! Real number of digits of output formats:
integer, parameter:: DR16P = 42   !< Number of digits of output format FR16P.
integer, parameter:: DR8P  = 23   !< Number of digits of output format FR8P.
integer, parameter:: DR4P  = 13   !< Number of digits of output format FR4P.
integer, parameter:: DR_P  = DR8P !< Number of digits of output format FR_P.
! Integer output formats:
character(5), parameter:: FI8P   = '(I20)'    !< Output format                     for kind=I8P variable.
character(8), parameter:: FI8PZP = '(I20.19)' !< Output format with zero prefixing for kind=I8P variable.
character(5), parameter:: FI4P   = '(I11)'    !< Output format                     for kind=I4P variable.
character(8), parameter:: FI4PZP = '(I11.10)' !< Output format with zero prefixing for kind=I4P variable.
character(4), parameter:: FI2P   = '(I6)'     !< Output format                     for kind=I2P variable.
character(6), parameter:: FI2PZP = '(I6.5)'   !< Output format with zero prefixing for kind=I2P variable.
character(4), parameter:: FI1P   = '(I4)'     !< Output format                     for kind=I1P variable.
character(6), parameter:: FI1PZP = '(I4.3)'   !< Output format with zero prefixing for kind=I1P variable.
character(5), parameter:: FI_P   = FI4P       !< Output format                     for kind=I_P variable.
character(8), parameter:: FI_PZP = FI4PZP     !< Output format with zero prefixing for kind=I_P variable.
! Integer number of digits of output formats:
integer, parameter:: DI8P = 20   !< Number of digits of output format I8P.
integer, parameter:: DI4P = 11   !< Number of digits of output format I4P.
integer, parameter:: DI2P = 6    !< Number of digits of output format I2P.
integer, parameter:: DI1P = 4    !< Number of digits of output format I1P.
integer, parameter:: DI_P = DI4P !< Number of digits of output format I_P.
! List of kinds
integer,       parameter:: NRknd=4                                           !< Number of defined real kinds.
integer,       parameter:: RPl(1:NRknd)=[R16P,R8P,R4P,R_P]                   !< List of defined real kinds.
character(10), parameter:: FRl(1:NRknd)=[FR16P,FR8P,FR4P//' ',FR_P]          !< List of defined real kinds output format.
integer,       parameter:: NIknd=5                                           !< Number of defined integer kinds.
integer,       parameter:: RIl(1:NIknd)=[I8P,I4P,I2P,I1P,I_P]                !< List of defined integer kinds.
character(5),  parameter:: FIl(1:NIknd)=[FI8P,FI4P,FI2P//' ',FI1P//' ',FI_P] !< List of defined integer kinds output format.

! Useful parameters for handling numbers ranges.
! Real min and max values:
real(R16P), parameter:: MinR16P = -huge(1._R16P), MaxR16P = huge(1._R16P) !< Min and max values of kind=R16P variable.
real(R8P),  parameter:: MinR8P  = -huge(1._R8P ), MaxR8P  = huge(1._R8P ) !< Min and max values of kind=R8P variable.
real(R4P),  parameter:: MinR4P  = -huge(1._R4P ), MaxR4P  = huge(1._R4P ) !< Min and max values of kind=R4P variable.
real(R_P),  parameter:: MinR_P  = MinR8P,         MaxR_P  = MaxR8P        !< Min and max values of kind=R_P variable.
! Real number of bits/bytes
integer(I2P):: BIR16P, BYR16P !< Number of bits/bytes of kind=R16P variable.
integer(I1P):: BIR8P,  BYR8P  !< Number of bits/bytes of kind=R8P variable.
integer(I1P):: BIR4P,  BYR4P  !< Number of bits/bytes of kind=R4P variable.
integer(I1P):: BIR_P,  BYR_P  !< Number of bits/bytes of kind=R_P variable.
! Real smallest values:
real(R16P), parameter:: smallR16P = tiny(1._R16P) !< Smallest representable value of kind=R16P variable.
real(R8P),  parameter:: smallR8P  = tiny(1._R8P ) !< Smallest representable value of kind=R8P variable.
real(R4P),  parameter:: smallR4P  = tiny(1._R4P ) !< Smallest representable value of kind=R4P variable.
real(R_P),  parameter:: smallR_P  = smallR8P      !< Smallest representable value of kind=R_P variable.
! Integer min and max values:
integer(I8P), parameter:: MinI8P = -huge(1_I8P), MaxI8P = huge(1_I8P) !< Min and max values of kind=I8P variable.
integer(I4P), parameter:: MinI4P = -huge(1_I4P), MaxI4P = huge(1_I4P) !< Min and max values of kind=I4P variable.
integer(I2P), parameter:: MinI2P = -huge(1_I2P), MaxI2P = huge(1_I2P) !< Min and max values of kind=I2P variable.
integer(I1P), parameter:: MinI1P = -huge(1_I1P), MaxI1P = huge(1_I1P) !< Min and max values of kind=I1P variable.
integer(I_P), parameter:: MinI_P = MinI4P,       MaxI_P = MaxI4P      !< Min and max values of kind=I_P variable.
! Integer number of bits/bytes:
integer(I8P), parameter:: BII8P = bit_size(MaxI8P), BYI8P = bit_size(MaxI8P)/8_I8P !< Number of bits/bytes of kind=I8P variable.
integer(I4P), parameter:: BII4P = bit_size(MaxI4P), BYI4P = bit_size(MaxI4P)/8_I4P !< Number of bits/bytes of kind=I4P variable.
integer(I2P), parameter:: BII2P = bit_size(MaxI2P), BYI2P = bit_size(MaxI2P)/8_I2P !< Number of bits/bytes of kind=I2P variable.
integer(I1P), parameter:: BII1P = bit_size(MaxI1P), BYI1P = bit_size(MaxI1P)/8_I1P !< Number of bits/bytes of kind=I1P variable.
integer(I_P), parameter:: BII_P = bit_size(MaxI_P), BYI_P = bit_size(MaxI_P)/8_I_P !< Number of bits/bytes of kind=I_P variable.
! Smallest real representable difference by the running calculator.
#ifdef pgf95
real(R16P), parameter:: ZeroR16 = 0._R16P
real(R8P),  parameter:: ZeroR8  = 0._R8P
real(R4P),  parameter:: ZeroR4  = 0._R4P
#else
real(R16P), parameter:: ZeroR16 = nearest(1._R16P, 1._R16P) - &
                                  nearest(1._R16P,-1._R16P) !< Smallest representable difference of kind=R16P variable.
real(R8P),  parameter:: ZeroR8  = nearest(1._R8P, 1._R8P) - &
                                  nearest(1._R8P,-1._R8P)   !< Smallest representable difference of kind=R8P variable.
real(R4P),  parameter:: ZeroR4  = nearest(1._R4P, 1._R4P) - &
                                  nearest(1._R4P,-1._R4P)   !< Smallest representable difference of kind=R4P variable.
#endif
real(R_P),  parameter:: Zero    = ZeroR8                    !< Smallest representable difference of kind=R_P variable.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Overloading of the intrinsic "bit_size" function for computing the number of bits of (also) real and character variables;
!> variable,     intent(\b IN)::  <b>\em n</b> input;
!> integer(I1P), intent(\b OUT):: <b>\em bits</b> output number of bits of input number.
!> @ingroup IR_PrecisionInterface
interface bit_size
  module procedure                &
#ifdef r16p
                   bit_size_R16p, &
#endif
                   bit_size_R8P,  &
                   bit_size_R4P,  &
                   bit_size_chr
endinterface
!> @brief Overloading of the "byte_size" function for computing the number of bytes.
!> @ingroup IR_PrecisionInterface
interface byte_size
  module procedure                 &
                   byte_size_I8P,  &
                   byte_size_I4P,  &
                   byte_size_I2P,  &
                   byte_size_I1P,  &
#ifdef r16p
                   byte_size_R16p, &
#endif
                   byte_size_R8P,  &
                   byte_size_R4P,  &
                   byte_size_chr
endinterface
!> @brief Procedure for converting number, real and integer, to string (number to string type casting);
!> logical, intent(\b IN), optional:: <b>\em no_sign</b> flag for do not write sign;
!> number,  intent(\b IN)::           <b>\em n</b> input number;
!> string,  intent(\b OUT)::          <b>\em str</b> output string.
!> @ingroup IR_PrecisionInterface
interface str
  module procedure                    &
#ifdef r16p
                   strf_R16P,str_R16P,&
#endif
                   strf_R8P ,str_R8P, &
                   strf_R4P ,str_R4P, &
                   strf_I8P ,str_I8P, &
                   strf_I4P ,str_I4P, &
                   strf_I2P ,str_I2P, &
                   strf_I1P ,str_I1P
endinterface
!> @brief Procedure for converting number, integer, to string, prefixing with the right number of zeros (number to string type
!> casting with zero padding);
!> number,  intent(\b IN), optional:: <b>\em no_zpad</b> number of padding zeros;
!> number,  intent(\b IN)::           <b>\em n  </b> input number;
!> string,  intent(\b OUT)::          <b>\em str</b> output string.
!> @ingroup IR_PrecisionInterface
interface strz
  module procedure strz_I8P,  &
                   strz_I4P,  &
                   strz_I2P,  &
                   strz_I1P
endinterface
!> @brief Procedure for converting string to number, real or initeger, (string to number type casting);
!> string,  intent(\b IN)::  <b>\em str</b> input string;
!> number,  intent(\b OUT):: <b>\em n  </b> output number.
!> @ingroup IR_PrecisionInterface
interface cton
  module procedure            &
#ifdef r16p
                   ctor_R16P, &
#endif
                   ctor_R8P,  &
                   ctor_R4P,  &
                   ctoi_I8P,  &
                   ctoi_I4P,  &
                   ctoi_I2P,  &
                   ctoi_I1P
endinterface
!> @brief Procedure for converting number, real and integer, to bit-string (number to bit-string type casting);
!> number,  intent(\b IN)::  <b>\em n</b> input number;
!> string,  intent(\b OUT):: <b>\em bstr</b> output bit-string.
!> @ingroup IR_PrecisionInterface
interface bstr
  module procedure           &
#ifdef r16p
                   bstr_R16P,&
#endif
                   bstr_R8P, &
                   bstr_R4P, &
                   bstr_I8P, &
                   bstr_I4P, &
                   bstr_I2P, &
                   bstr_I1P
endinterface
!> @brief Procedure for converting bit-string to number, real or initeger, (bit-string to number type casting);
!> string,  intent(\b IN)::  <b>\em bstr</b> input bit-string;
!> number,  intent(\b OUT):: <b>\em n</b> output number.
!> @ingroup IR_PrecisionInterface
interface bcton
  module procedure            &
#ifdef r16p
                   bctor_R16P, &
#endif
                   bctor_R8P,  &
                   bctor_R4P,  &
                   bctoi_I8P,  &
                   bctoi_I4P,  &
                   bctoi_I2P,  &
                   bctoi_I1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup IR_PrecisionPublicProcedure
  !> @{
  !> @brief Procedure for checking if the type of the bit ordering of the running architecture is little endian.
  pure function is_little_endian() result(is_little)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical::      is_little !< Logical output: true is the running architecture uses little endian ordering, false otherwise.
  integer(I1P):: int1(1:4) !< One byte integer array for casting 4 bytes integer.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  int1 = transfer(1_I4P,int1)
  is_little = (int1(1)==1_I1P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_little_endian

  !> @brief Subroutine for checking the type of bit ordering (big or little endian) of the running architecture; the result is
  !> stored into the "endian" global variable.
  subroutine check_endian()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (is_little_endian()) then
    endian = endianL
  else
    endian = endianB
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine check_endian
  !> @}

  !> @ingroup IR_PrecisionPrivateProcedure
  !> @{
  !> @brief Procedure for computing the number of bits of a real variable.
  elemental function bit_size_R16P(r) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: r       !< Real variable whose number of bits must be computed.
  integer(I2P)::           bits    !< Number of bits of r.
  integer(I1P)::           mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(r,mold),dim=1,kind=I2P)*8_I2P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_R16P

  !> @brief Procedure for computing the number of bits of a real variable.
  elemental function bit_size_R8P(r) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: r       !< Real variable whose number of bits must be computed.
  integer(I1P)::          bits    !< Number of bits of r.
  integer(I1P)::          mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(r,mold),dim=1,kind=I1P)*8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_R8P

  !> @brief Procedure for computing the number of bits of a real variable.
  elemental function bit_size_R4P(r) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: r       !< Real variable whose number of bits must be computed.
  integer(I1P)::          bits    !< Number of bits of r.
  integer(I1P)::          mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(r,mold),dim=1,kind=I1P)*8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_R4P

  !> @brief Procedure for computing the number of bits of a character variable.
  elemental function bit_size_chr(c) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: c       !< Character variable whose number of bits must be computed.
  integer(I4P)::             bits    !< Number of bits of c.
  integer(I1P)::             mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(c,mold),dim=1,kind=I4P)*8_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_chr

  !> @brief Procedure for computing the number of bytes of an integer variable.
  elemental function byte_size_I8P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)::             bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I8P

  !> @brief Procedure for computing the number of bytes of an integer variable.
  elemental function byte_size_I4P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)::             bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I4P

  !> @brief Procedure for computing the number of bytes of an integer variable.
  elemental function byte_size_I2P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN):: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)::             bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I2P

  !> @brief Procedure for computing the number of bytes of an integer variable.
  elemental function byte_size_I1P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)::             bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I1P

  !> @brief Procedure for computing the number of bytes of a real variable.
  elemental function byte_size_R16P(r) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: r     !< Real variable whose number of bytes must be computed.
  integer(I1P)::           bytes !< Number of bytes of r.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(r)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_R16P

  !> @brief Procedure for computing the number of bytes of a real variable.
  elemental function byte_size_R8P(r) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: r     !< Real variable whose number of bytes must be computed.
  integer(I1P)::          bytes !< Number of bytes of r.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(r)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_R8P

  !> @brief Procedure for computing the number of bytes of a real variable.
  elemental function byte_size_R4P(r) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: r     !< Real variable whose number of bytes must be computed.
  integer(I1P)::          bytes !< Number of bytes of r.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(r)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_R4P

  !> @brief Procedure for computing the number of bytes of a character variable.
  elemental function byte_size_chr(c) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: c     !< Character variable whose number of bytes must be computed.
  integer(I4P)::             bytes !< Number of bytes of c.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(c)/8_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_chr

  !> @brief Procedure for converting real to string. This function achieves casting of real to string.
  elemental function strf_R16P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  real(R16P),   intent(IN):: n   !< Real to be converted.
  character(DR16P)::         str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_R16P

  !> @brief Procedure for converting real to string. This function achieves casting of real to string.
  elemental function strf_R8P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  real(R8P),    intent(IN):: n   !< Real to be converted.
  character(DR8P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_R8P

  !> @brief Procedure for converting real to string. This function achieves casting of real to string.
  elemental function strf_R4P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  real(R4P),    intent(IN):: n   !< Real to be converted.
  character(DR4P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_R4P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function strf_I8P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  integer(I8P), intent(IN):: n   !< Integer to be converted.
  character(DI8P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I8P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function strf_I4P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  integer(I4P), intent(IN):: n   !< Integer to be converted.
  character(DI4P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I4P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function strf_I2P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  integer(I2P), intent(IN):: n   !< Integer to be converted.
  character(DI2P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I2P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function strf_I1P(fm,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: fm  !< Format different from the standard for the kind.
  integer(I1P), intent(IN):: n   !< Integer to be converted.
  character(DI1P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,trim(fm)) n ! Casting of n to string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I1P

  !> @brief Procedure for converting real to string. This function achieves casting of real to string.
  elemental function str_R16P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,    intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  real(R16P), intent(IN)::           n       !< Real to be converted.
  character(DR16P)::                 str     !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FR16P) n                ! Casting of n to string.
  if (n>0._R16P) str(1:1)='+'       ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R16P

  !> @brief Procedure for converting real to string. This function achieves casting of real to string.
  elemental function str_R8P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,      intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  real(R8P),    intent(IN)::           n       !< Real to be converted.
  character(DR8P)::                    str     !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FR8P) n                 ! Casting of n to string.
  if (n>0._R8P) str(1:1)='+'        ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R8P

  !> @brief Procedure for converting real to string. This function achieves casting of real to string.
  elemental function str_R4P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,   intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  real(R4P), intent(IN)::           n       !< Real to be converted.
  character(DR4P)::                 str     !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FR4P) n                 ! Casting of n to string.
  if (n>0._R4P) str(1:1)='+'        ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R4P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function str_I8P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,      intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  integer(I8P), intent(IN)::           n       !< Integer to be converted.
  character(DI8P)::                    str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI8P) n                 ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I8P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I8P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function str_I4P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,      intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  integer(I4P), intent(IN)::           n       !< Integer to be converted.
  character(DI4P)::                    str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI4P) n                 ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I4P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I4P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function str_I2P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,      intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  integer(I2P), intent(IN)::           n       !< Integer to be converted.
  character(DI2P)::                    str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI2P) n                 ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I2P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I2P

  !> @brief Procedure for converting integer to string. This function achieves casting of integer to string.
  elemental function str_I1P(no_sign,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,      intent(IN), optional:: no_sign !< Flag for leaving out the sign.
  integer(I1P), intent(IN)::           n       !< Integer to be converted.
  character(DI1P)::                    str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI1P) n                 ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I1P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I1P

  !> @brief Procedure for converting integer to string, prefixing with the right number of zeros. This function achieves casting of
  !> integer to string.
  elemental function strz_I8P(nz_pad,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: nz_pad !< Number of zeros padding.
  integer(I8P), intent(IN)::           n      !< Integer to be converted.
  character(DI8P)::                    str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI8PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI8P-nz_pad:DI8P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I8P

  !> @brief Procedure for converting integer to string, prefixing with the right number of zeros. This function achieves casting of
  !> integer to string.
  elemental function strz_I4P(nz_pad,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: nz_pad !< Number of zeros padding.
  integer(I4P), intent(IN)::           n      !< Integer to be converted.
  character(DI4P)::                    str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI4PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI4P-nz_pad:DI4P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I4P

  !> @brief Procedure for converting integer to string, prefixing with the right number of zeros. This function achieves casting of
  !> integer to string.
  elemental function strz_I2P(nz_pad,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: nz_pad !< Number of zeros padding.
  integer(I2P), intent(IN)::           n      !< Integer to be converted.
  character(DI2P)::                    str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI2PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI2P-nz_pad:DI2P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I2P

  !> @brief Procedure for converting integer to string, prefixing with the right number of zeros. This function achieves casting of
  !> integer to string.
  elemental function strz_I1P(nz_pad,n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: nz_pad !< Number of zeros padding.
  integer(I1P), intent(IN)::           n      !< Integer to be converted.
  character(DI1P)::                    str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI1PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI1P-nz_pad:DI1P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I1P

  !> @brief Procedure for converting string to real. This function achieves casting of string to real.
  function ctor_R16P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  real(R16P),   intent(IN):: knd !< Number kind.
  real(R16P)::               n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')             'Conversion of string "'//str//'" to real failed'
    write(stderr,'(A,'//FR16P//')') 'Kind parameter ',knd
    write(stderr,'(A)')             'Function used "ctor_R16P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctor_R16P

  !> @brief Procedure for converting string to real. This function achieves casting of string to real.
  function ctor_R8P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  real(R8P),    intent(IN):: knd !< Number kind.
  real(R8P)::                n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')            'Conversion of string "'//str//'" to real failed'
    write(stderr,'(A,'//FR8P//')') 'Kind parameter ',knd
    write(stderr,'(A)')            'Function used "ctor_R8P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctor_R8P

  !> @brief Procedure for converting string to real. This function achieves casting of string to real.
  function ctor_R4P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  real(R4P),    intent(IN):: knd !< Number kind.
  real(R4P)::                n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')            'Conversion of string "'//str//'" to real failed'
    write(stderr,'(A,'//FR4P//')') 'Kind parameter ',knd
    write(stderr,'(A)')            'Function used "ctor_R4P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctor_R4P

  !> @brief Procedure for converting string to integer. This function achieves casting of string to integer.
  function ctoi_I8P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  integer(I8P), intent(IN):: knd !< Number kind.
  integer(I8P)::             n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')            'Conversion of string "'//str//'" to integer failed'
    write(stderr,'(A,'//FI8P//')') 'Kind parameter ',knd
    write(stderr,'(A)')            'Function used "ctoi_I8P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I8P

  !> @brief Procedure for converting string to integer. This function achieves casting of string to integer.
  function ctoi_I4P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  integer(I4P), intent(IN):: knd !< Number kind.
  integer(I4P)::             n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')            'Conversion of string "'//str//'" to integer failed'
    write(stderr,'(A,'//FI4P//')') 'Kind parameter ',knd
    write(stderr,'(A)')            'Function used "ctoi_I4P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I4P

  !> @brief Procedure for converting string to integer. This function achieves casting of string to integer.
  function ctoi_I2P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  integer(I2P), intent(IN):: knd !< Number kind.
  integer(I2P)::             n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')            'Conversion of string "'//str//'" to integer failed'
    write(stderr,'(A,'//FI2P//')') 'Kind parameter ',knd
    write(stderr,'(A)')            'Function used "ctoi_I2P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I2P

  !> @brief Procedure for converting string to integer. This function achieves casting of string to integer.
  function ctoi_I1P(str,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: str !< String containing input number.
  integer(I1P), intent(IN):: knd !< Number kind.
  integer(I1P)::             n   !< Number returned.
  integer(I4P)::             err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str,*,iostat=err) n ! Casting of str to n.
  if (err/=0) then
    write(stderr,'(A)')            'Conversion of string "'//str//'" to integer failed'
    write(stderr,'(A,'//FI1P//')') 'Kind parameter ',knd
    write(stderr,'(A)')            'Function used "ctoi_I1P"'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I1P

  !> @brief Procedure for converting real to string of bits. This function achieves casting of real to bit-string.
  !> @note It is assumed that R16P is represented by means of 128 bits, but this is not ensured in all architectures.
  elemental function bstr_R16P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: n    !< Real to be converted.
  character(128)::        bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B128.128)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_R16P

  !> @brief Procedure for converting real to string of bits. This function achieves casting of real to bit-string.
  !> @note It is assumed that R8P is represented by means of 64 bits, but this is not ensured in all architectures.
  elemental function bstr_R8P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: n    !< Real to be converted.
  character(64)::         bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B64.64)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_R8P

  !> @brief Procedure for converting real to string of bits. This function achieves casting of real to bit-string.
  !> @note It is assumed that R4P is represented by means of 32 bits, but this is not ensured in all architectures.
  elemental function bstr_R4P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: n    !< Real to be converted.
  character(32)::         bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B32.32)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_R4P

  !> @brief Procedure for converting integer to string of bits. This function achieves casting of integer to bit-string.
  !> @note It is assumed that I8P is represented by means of 64 bits, but this is not ensured in all architectures.
  elemental function bstr_I8P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: n    !< Real to be converted.
  character(64)::            bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B64.64)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I8P

  !> @brief Procedure for converting integer to string of bits. This function achieves casting of integer to bit-string.
  !> @note It is assumed that I4P is represented by means of 32 bits, but this is not ensured in all architectures.
  elemental function bstr_I4P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: n    !< Real to be converted.
  character(32)::            bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B32.32)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I4P

  !> @brief Procedure for converting integer to string of bits. This function achieves casting of integer to bit-string.
  !> @note It is assumed that I2P is represented by means of 16 bits, but this is not ensured in all architectures.
  elemental function bstr_I2P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN):: n    !< Real to be converted.
  character(16)::            bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B16.16)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I2P

  !> @brief Procedure for converting integer to string of bits. This function achieves casting of integer to bit-string.
  !> @note It is assumed that I1P is represented by means of 8 bits, but this is not ensured in all architectures.
  elemental function bstr_I1P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: n    !< Real to be converted.
  character(8)::             bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr,'(B8.8)')n ! Casting of n to bit-string.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I1P

  !> @brief Procedure for converting bit-string to real. This function achieves casting of bit-string to real.
  elemental function bctor_R8P(bstr,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: bstr !< String containing input number.
  real(R8P),    intent(IN):: knd  !< Number kind.
  real(R8P)::                n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(.true.,bit_size(knd)))//'.'//trim(str(.true.,bit_size(knd)))//')')n ! Casting of bstr to n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctor_R8P

  !> @brief Procedure for converting bit-string to real. This function achieves casting of bit-string to real.
  elemental function bctor_R4P(bstr,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: bstr !< String containing input number.
  real(R4P),    intent(IN):: knd  !< Number kind.
  real(R4P)::                n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(.true.,bit_size(knd)))//'.'//trim(str(.true.,bit_size(knd)))//')')n ! Casting of bstr to n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctor_R4P

  !> @brief Procedure for converting bit-string to integer. This function achieves casting of bit-string to integer.
  elemental function bctoi_I8P(bstr,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: bstr !< String containing input number.
  integer(I8P), intent(IN):: knd  !< Number kind.
  integer(I8P)::             n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(.true.,bit_size(knd)))//'.'//trim(str(.true.,bit_size(knd)))//')')n ! Casting of bstr to n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I8P

  !> @brief Procedure for converting bit-string to integer. This function achieves casting of bit-string to integer.
  elemental function bctoi_I4P(bstr,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: bstr !< String containing input number.
  integer(I4P), intent(IN):: knd  !< Number kind.
  integer(I4P)::             n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(.true.,bit_size(knd)))//'.'//trim(str(.true.,bit_size(knd)))//')')n ! Casting of bstr to n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I4P

  !> @brief Procedure for converting bit-string to integer. This function achieves casting of bit-string to integer.
  elemental function bctoi_I2P(bstr,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: bstr !< String containing input number.
  integer(I2P), intent(IN):: knd  !< Number kind.
  integer(I2P)::             n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(.true.,bit_size(knd)))//'.'//trim(str(.true.,bit_size(knd)))//')')n ! Casting of bstr to n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I2P

  !> @brief Procedure for converting bit-string to integer. This function achieves casting of bit-string to integer.
  elemental function bctoi_I1P(bstr,knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: bstr !< String containing input number.
  integer(I1P), intent(IN):: knd  !< Number kind.
  integer(I1P)::             n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(.true.,bit_size(knd)))//'.'//trim(str(.true.,bit_size(knd)))//')')n ! Casting of bstr to n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I1P
  !> @}

  !> Subroutine for initilizing module's variables that are not initialized into the definition specification.
  !> @ingroup IR_PrecisionPublicProcedure
  subroutine IR_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! checking the bit ordering architecture
  call check_endian
  ! computing the bits/bytes sizes of real variables
  BIR8P  = bit_size(r=MaxR8P)  ; BYR8P  = BIR8P/8_I1P
  BIR4P  = bit_size(r=MaxR4P)  ; BYR4P  = BIR4P/8_I1P
  BIR_P  = bit_size(r=MaxR_P)  ; BYR_P  = BIR_P/8_I1P
#ifdef r16p
  BIR16P = bit_size(r=MaxR16P) ; BYR16P = BIR16P/8_I2P
#else
  BIR16P = int(BIR8P,kind=I2P) ; BYR16P = BIR16P/8_I2P
#endif
  ir_initialized = .true.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IR_init

  !>Subroutine for printing to the standard output the kind definition of reals and integers and the utility variables.
  !> @ingroup IR_PrecisionPublicProcedure
  subroutine IR_Print(myrank,Nproc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: myrank  !< Actual rank process necessary for concurrent multi-processes calls.
  integer(I4P), intent(IN), optional:: Nproc   !< Number of MPI processes used.
  character(DI4P)::                    rks     !< String containing myrank.
  integer(I4P)::                       rank,Np !< Dummy temporary variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.ir_initialized) call IR_init
  rank = 0 ; if (present(myrank)) rank = myrank ; Np = 1 ; if (present(Nproc)) Np = Nproc ; rks = 'rank'//trim(strz(Np,rank))
  ! printing informations
  if (endian==endianL) then
    write(stdout,'(A)')        trim(rks)//' This architecture has LITTLE Endian bit ordering'
  else
    write(stdout,'(A)')        trim(rks)//' This architecture has BIG Endian bit ordering'
  endif
  write(stdout,'(A)')          trim(rks)//' Reals kind precision definition'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R16P Kind "',R16P,'" | FR16P format "'//FR16P//'" | DR16P chars "',DR16P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R8P  Kind "',R8P, '" | FR8P  format "'//FR8P// '" | DR8P  chars "',DR8P ,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R4P  Kind "',R4P, '" | FR4P  format "'//FR4P//'"  | DR4P  chars "',DR4P ,'"'
  write(stdout,'(A)')          trim(rks)//' Integers kind precision definition'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I8P Kind "',I8P,'" | FI8P format "'//FI8P// '" | DI8P chars "',DI8P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I4P Kind "',I4P,'" | FI4P format "'//FI4P// '" | DI4P chars "',DI4P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I2P Kind "',I2P,'" | FI2P format "'//FI2P//'"  | DI2P chars "',DI2P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I1P Kind "',I1P,'" | FI1P format "'//FI1P//'"  | DI1P chars "',DI1P,'"'
  write(stdout,'(A)')          trim(rks)//' Reals minimum and maximum values'
  write(stdout,'(A)')          trim(rks)//' MinR16P "'//trim(str(n=MinR16P))//'" | MaxR16P "'//trim(str(n=MaxR16P))//'"'
  write(stdout,'(A)')          trim(rks)//' MinR8P  "'//trim(str(n=MinR8P))// '" | MaxR8P  "'//trim(str(n=MaxR8P))// '"'
  write(stdout,'(A)')          trim(rks)//' MinR4P  "'//trim(str(n=MinR4P))// '" | MaxR4P  "'//trim(str(n=MaxR4P))// '"'
  write(stdout,'(A)')          trim(rks)//' Reals bits/bytes sizes'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R16P bits "',BIR16P,'", bytes "',BYR16P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R8P bits "', BIR8P, '", bytes "',BYR8P, '"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R4P bits "', BIR4P, '", bytes "',BYR4P, '"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' R_P bits "', BIR_P, '", bytes "',BYR_P, '"'
  write(stdout,'(A)')          trim(rks)//' Integers minimum and maximum values'
  write(stdout,'(A)')          trim(rks)//' MinI8P "'//trim(str(n=MinI8P))//'" | MaxI8P "'//trim(str(n=MaxI8P))//'"'
  write(stdout,'(A)')          trim(rks)//' MinI4P "'//trim(str(n=MinI4P))//'" | MaxI4P "'//trim(str(n=MaxI4P))//'"'
  write(stdout,'(A)')          trim(rks)//' MinI2P "'//trim(str(n=MinI2P))//'" | MaxI2P "'//trim(str(n=MaxI2P))//'"'
  write(stdout,'(A)')          trim(rks)//' MinI1P "'//trim(str(n=MinI1P))//'" | MaxI1P "'//trim(str(n=MaxI1P))//'"'
  write(stdout,'(A)')          trim(rks)//' Integers bits/bytes sizes'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I8P bits "',BII8P,'", bytes "',BYI8P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I4P bits "',BII4P,'", bytes "',BYI4P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I2P bits "',BII2P,'", bytes "',BYI2P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I1P bits "',BII1P,'", bytes "',BYI1P,'"'
  write(stdout,'(A,I2,A,I2,A)')trim(rks)//' I_P bits "',BII_P,'", bytes "',BYI_P,'"'
  write(stdout,'(A)')          trim(rks)//' Machine precisions'
  write(stdout,'(A)')          trim(rks)//' ZeroR16 "'//trim(str(.true.,ZeroR16))//'"'
  write(stdout,'(A)')          trim(rks)//' ZeroR8  "'//trim(str(.true.,ZeroR8 ))//'"'
  write(stdout,'(A)')          trim(rks)//' ZeroR4  "'//trim(str(.true.,ZeroR4 ))//'"'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IR_Print
endmodule IR_Precision
