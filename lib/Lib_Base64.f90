!> @ingroup Library
!> @{
!> @defgroup Lib_Base64Library Lib_Base64
!> base64 encoding/decoding library
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Lib_Base64Interface Lib_Base64
!> base64 encoding/decoding library
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Base64PublicProcedure Lib_Base64
!> base64 encoding/decoding library
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Base64PrivateProcedure Lib_Base64
!> base64 encoding/decoding library
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_Base64GlobalVarPar Lib_Base64
!> base64 encoding/decoding library
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_Base64PrivateVarPar Lib_Base64
!> base64 encoding/decoding library
!> @}

!> This module contains base64 encoding/decoding procedures.
!> @todo \b Decoding: Implement decoding functions.
!> @todo \b DocComplete: Complete the documentation.
!> @ingroup Lib_Base64Library
module Lib_Base64
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: b64_encode
!public:: b64_decode
public:: pack_data
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Lib_Base64GlobalVarPar
!> @{
!> @}
!> @ingroup Lib_Base64PrivateVarPar
!> @{
character(64):: base64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/" !< Base64 alphabet.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Subroutine for encoding numbers (integer and real) to base64.
!> @ingroup Lib_Base64Interface
interface b64_encode
  module procedure b64_encode_R8_a, &
                   b64_encode_R4_a, &
                   b64_encode_I8_a, &
                   b64_encode_I4_a, &
                   b64_encode_I2_a, &
                   b64_encode_I1_a
endinterface
!!> @brief Subroutine for decoding numbers (integer and real) from base64.
!!> @ingroup Lib_Base64Interface
!interface b64_decode
!  module procedure b64_decode_R8_a, &
!                   b64_decode_R4_a, &
!                   b64_decode_I8_a, &
!                   b64_decode_I4_a, &
!                   b64_decode_I2_a, &
!                   b64_decode_I1_a
!endinterface
!> @brief Subroutine for packing different kinds of data into single I1P array. This is useful for encoding different kinds
!> variables into a single stream of bits.
!> @ingroup Lib_Base64Interface
interface pack_data
  module procedure pack_data_R8_R4,pack_data_R8_I8,pack_data_R8_I4,pack_data_R8_I2,pack_data_R8_I1, &
                   pack_data_R4_R8,pack_data_R4_I8,pack_data_R4_I4,pack_data_R4_I2,pack_data_R4_I1, &
                   pack_data_I8_R8,pack_data_I8_R4,pack_data_I8_I4,pack_data_I8_I2,pack_data_I8_I1, &
                   pack_data_I4_R8,pack_data_I4_R4,pack_data_I4_I8,pack_data_I4_I2,pack_data_I4_I1, &
                   pack_data_I2_R8,pack_data_I2_R4,pack_data_I2_I8,pack_data_I2_I4,pack_data_I2_I1, &
                   pack_data_I1_R8,pack_data_I1_R4,pack_data_I1_I8,pack_data_I1_I4,pack_data_I1_I2
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Base64PrivateProcedure
  !> @{
  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R8_R4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),                 intent(IN)::    a1(1:)    !< Firs data stream.
  real(R4P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R8_R4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R8_I8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I8P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R8_I8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R8_I4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I4P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R8_I4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R8_I2(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I2P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R8_I2

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R8_I1(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I1P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R8_I1

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R4_R8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),                 intent(IN)::    a1(1:)    !< Firs data stream.
  real(R8P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R4_R8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R4_I8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I8P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R4_I8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R4_I4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I4P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R4_I4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R4_I2(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I2P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R4_I2

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_R4_I1(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),                 intent(IN)::    a1(1:)    !< First data stream.
  integer(I1P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_R4_I1

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I8_R8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R8P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I8_R8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I8_R4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R4P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I8_R4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I8_I4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I4P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I8_I4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I8_I2(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I2P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I8_I2

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I8_I1(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I1P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I8_I1

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I4_R8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R8P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I4_R8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I4_R4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R4P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I4_R4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I4_I8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I8P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I4_I8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I4_I2(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I2P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I4_I2

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I4_I1(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I1P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I4_I1

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I2_R8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R8P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I2_R8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I2_R4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R4P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I2_R4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I2_I8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I8P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I2_I8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I2_I4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I4P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I2_I4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I2_I1(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I1P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I2_I1

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I1_R8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R8P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I1_R8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I1_R4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),              intent(IN)::    a1(1:)    !< First data stream.
  real(R4P),                 intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I1_R4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I1_I8(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I8P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I1_I8

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I1_I4(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I4P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I1_I4

  !> @brief Subroutine for packing different kinds of data into single I1P array.
  pure subroutine pack_data_I1_I2(a1,a2,packed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),              intent(IN)::    a1(1:)    !< First data stream.
  integer(I2P),              intent(IN)::    a2(1:)    !< Second data stream.
  integer(I1P), allocatable, intent(INOUT):: packed(:) !< Packed data into I1P array.
  integer(I1P), allocatable::                p1(:)     !< Temporary packed data of first stream.
  integer(I1P), allocatable::                p2(:)     !< Temporary packed data of second stream.
  integer(I4P)::                             np        !< Size of temporary packed data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  np = size(transfer(a1,p1)) ; allocate(p1(1:np)) ; p1 = transfer(a1,p1)
  np = size(transfer(a2,p2)) ; allocate(p2(1:np)) ; p2 = transfer(a2,p2)
  if (allocated(packed)) deallocate(packed) ; allocate(packed(1:size(p1,dim=1)+size(p2,dim=1))) ; packed = [p1,p2]
  deallocate(p1,p2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pack_data_I1_I2

  !> @brief Subroutine for encoding bits (must be multiple of 24 bits) into base64 charcaters code (of length multiple of 4).
  !> @note The bits stream are encoded in chunks of 24 bits as the following example (in little endian order):
  !> @code
  !> +--first octet--+-second octet--+--third octet--+
  !> |7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|
  !> +-----------+---+-------+-------+---+-----------+
  !> |5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|
  !> +--1.index--+--2.index--+--3.index--+--4.index--+
  !> @endcode
  !> The 4 indexes are stored into 4 elements 8 bits array, thus 2 bits of each array element are not used.
  pure subroutine encode_bits(bits,padd,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN)::  bits(1:)  !< Bits to be encoded.
  integer(I4P), intent(IN)::  padd      !< Number of padding characters ('=').
  character(1), intent(OUT):: code(1:)  !< Characters code.
  integer(I1P)::              sixb(1:4) !< 6 bits slices (stored into 8 bits integer) of 24 bits input.
  integer(I8P)::              c,e       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c = 1_I8P
  do e=1_I8P,size(bits,dim=1),3_I8P ! loop over array elements: 3 bytes (24 bits) scanning
    sixb = 0_I1P
    call mvbits(bits(e  ),2,6,sixb(1),0)
    call mvbits(bits(e  ),0,2,sixb(2),4) ; call mvbits(bits(e+1),4,4,sixb(2),0)
    call mvbits(bits(e+1),0,4,sixb(3),2) ; call mvbits(bits(e+2),6,2,sixb(3),0)
    call mvbits(bits(e+2),0,6,sixb(4),0)
    sixb = sixb + 1_I1P
    code(c  :c  )(1:1) = base64(sixb(1):sixb(1))
    code(c+1:c+1)(1:1) = base64(sixb(2):sixb(2))
    code(c+2:c+2)(1:1) = base64(sixb(3):sixb(3))
    code(c+3:c+3)(1:1) = base64(sixb(4):sixb(4))
    c = c + 4_I8P
  enddo
  if (padd>0) code(size(code,dim=1)-padd+1:)(1:1)='='
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine encode_bits

  !> @brief Subroutine for encoding array numbers to base64 (R8P).
  pure subroutine b64_encode_R8_a(nB,n,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  nB                                 !< Number of bytes of single element of n.
  real(R8P),                 intent(IN)::  n(1:)                              !< Array of numbers to be encoded.
  character(1), allocatable, intent(OUT):: code(:)                            !< Encoded array.
  integer(I1P)::                           nI1P(1:((size(n,dim=1)*nB+2)/3)*3) !< One byte integer array containing n.
  integer(I4P)::                           padd                               !< Number of padding characters ('=').
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(code)) deallocate(code) ; allocate(code(1:((size(n,dim=1)*nB+2)/3)*4)) ! allocating code chars
  nI1P = transfer(n,nI1P)                                                              ! casting n to integer array of 1 byte elem
  padd = mod((size(n,dim=1)*nB),3_I4P) ; if (padd>0_I4P) padd = 3_I4P - padd           ! computing the number of padding characters
  call encode_bits(bits=nI1P,padd=padd,code=code)                                      ! encoding bits
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine b64_encode_R8_a

  !> @brief Subroutine for encoding array numbers to base64 (R4P).
  pure subroutine b64_encode_R4_a(nB,n,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  nB                                 !< Number of bytes of single element of n.
  real(R4P),                 intent(IN)::  n(1:)                              !< Array of numbers to be encoded.
  character(1), allocatable, intent(OUT):: code(:)                            !< Encoded array.
  integer(I1P)::                           nI1P(1:((size(n,dim=1)*nB+2)/3)*3) !< One byte integer array containing n.
  integer(I4P)::                           padd                               !< Number of padding characters ('=').
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(code)) deallocate(code) ; allocate(code(1:((size(n,dim=1)*nB+2)/3)*4)) ! allocating code chars
  nI1P = transfer(n,nI1P)                                                              ! casting n to integer array of 1 byte elem
  padd = mod((size(n,dim=1)*nB),3_I4P) ; if (padd>0_I4P) padd = 3_I4P - padd           ! computing the number of padding characters
  call encode_bits(bits=nI1P,padd=padd,code=code)                                      ! encoding bits
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine b64_encode_R4_a

  !> @brief Subroutine for encoding array numbers to base64 (I8P).
  pure subroutine b64_encode_I8_a(nB,n,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  nB                                 !< Number of bytes of single element of n.
  integer(I8P),              intent(IN)::  n(1:)                              !< Array of numbers to be encoded.
  character(1), allocatable, intent(OUT):: code(:)                            !< Encoded array.
  integer(I1P)::                           nI1P(1:((size(n,dim=1)*nB+2)/3)*3) !< One byte integer array containing n.
  integer(I4P)::                           padd                               !< Number of padding characters ('=').
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(code)) deallocate(code) ; allocate(code(1:((size(n,dim=1)*nB+2)/3)*4)) ! allocating code chars
  nI1P = transfer(n,nI1P)                                                              ! casting n to integer array of 1 byte elem
  padd = mod((size(n,dim=1)*nB),3_I4P) ; if (padd>0_I4P) padd = 3_I4P - padd           ! computing the number of padding characters
  call encode_bits(bits=nI1P,padd=padd,code=code)                                      ! encoding bits
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine b64_encode_I8_a

  !> @brief Subroutine for encoding array numbers to base64 (I4P).
  pure subroutine b64_encode_I4_a(nB,n,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  nB                                 !< Number of bytes of single element of n.
  integer(I4P),              intent(IN)::  n(1:)                              !< Array of numbers to be encoded.
  character(1), allocatable, intent(OUT):: code(:)                            !< Encoded array.
  integer(I1P)::                           nI1P(1:((size(n,dim=1)*nB+2)/3)*3) !< One byte integer array containing n.
  integer(I4P)::                           padd                               !< Number of padding characters ('=').
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(code)) deallocate(code) ; allocate(code(1:((size(n,dim=1)*nB+2)/3)*4)) ! allocating code chars
  nI1P = transfer(n,nI1P)                                                              ! casting n to integer array of 1 byte elem
  padd = mod((size(n,dim=1)*nB),3_I4P) ; if (padd>0_I4P) padd = 3_I4P - padd           ! computing the number of padding characters
  call encode_bits(bits=nI1P,padd=padd,code=code)                                      ! encoding bits
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine b64_encode_I4_a

  !> @brief Subroutine for encoding array numbers to base64 (I2P).
  pure subroutine b64_encode_I2_a(nB,n,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  nB                                 !< Number of bytes of single element of n.
  integer(I2P),              intent(IN)::  n(1:)                              !< Array of numbers to be encoded.
  character(1), allocatable, intent(OUT):: code(:)                            !< Encoded array.
  integer(I1P)::                           nI1P(1:((size(n,dim=1)*nB+2)/3)*3) !< One byte integer array containing n.
  integer(I4P)::                           padd                               !< Number of padding characters ('=').
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(code)) deallocate(code) ; allocate(code(1:((size(n,dim=1)*nB+2)/3)*4)) ! allocating code chars
  nI1P = transfer(n,nI1P)                                                              ! casting n to integer array of 1 byte elem
  padd = mod((size(n,dim=1)*nB),3_I4P) ; if (padd>0_I4P) padd = 3_I4P - padd           ! computing the number of padding characters
  call encode_bits(bits=nI1P,padd=padd,code=code)                                      ! encoding bits
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine b64_encode_I2_a

  !> @brief Subroutine for encoding array numbers to base64 (I1P).
  pure subroutine b64_encode_I1_a(nB,n,code)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  nB                                 !< Number of bytes of single element of n.
  integer(I1P),              intent(IN)::  n(1:)                              !< Array of numbers to be encoded.
  character(1), allocatable, intent(OUT):: code(:)                            !< Encoded array.
  integer(I1P)::                           nI1P(1:((size(n,dim=1)*nB+2)/3)*3) !< One byte integer array containing n.
  integer(I4P)::                           padd                               !< Number of padding characters ('=').
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(code)) deallocate(code) ; allocate(code(1:((size(n,dim=1)*nB+2)/3)*4)) ! allocating code chars
  nI1P = transfer(n,nI1P)                                                              ! casting n to integer array of 1 byte elem
  padd = mod((size(n,dim=1)*nB),3_I4P) ; if (padd>0_I4P) padd = 3_I4P - padd           ! computing the number of padding characters
  call encode_bits(bits=nI1P,padd=padd,code=code)                                      ! encoding bits
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine b64_encode_I1_a

  !!> @brief Subroutine for decoding array numbers from base64 (R8P).
  !pure subroutine b64_decode_R8_a(code,n)
  !!--------------------------------------------------------------------------------------------------------------------------------
  !implicit none
  !real(R8P),                      intent(OUT):: n(1:) !< Number to be decoded.
  !character(ncR8P*size(n,dim=1)), intent(IN)::  code  !< Encoded number.
  !integer(I4P)::                                c,d   !< Counters.
  !!--------------------------------------------------------------------------------------------------------------------------------

  !!--------------------------------------------------------------------------------------------------------------------------------
  !d = 1_I4P
  !do c=1,len(code),ncR8P
  !  call b64_decode_R8_s(code=code(c:c+ncR8P-1),n=n(d))
  !  d = d + 1_I4P
  !enddo
  !return
  !!--------------------------------------------------------------------------------------------------------------------------------
  !endsubroutine b64_decode_R8_a
  !> @}
endmodule Lib_Base64
