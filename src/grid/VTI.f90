!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Read data from image files of the visualization toolkit.
!--------------------------------------------------------------------------------------------------
module VTI
  use prec
  use zlib
  use base64
  use IO

  implicit none(type,external)
  private

  public :: &
    VTI_readDataset_real, &
    VTI_readDataset_int, &
    VTI_readGeometry

contains

!--------------------------------------------------------------------------------------------------
!> @brief Read integer dataset from a VTK image data (*.vti) file.
!> @details https://vtk.org/Wiki/VTK_XML_Formats
!--------------------------------------------------------------------------------------------------
function VTI_readDataset_int(fileContent,label) result(dataset)

  character(len=*), intent(in) :: &
    label, &
    fileContent
  integer, dimension(:), allocatable :: &
    dataset

  character(len=:), allocatable :: dataType, headerType, base64Str
  logical :: compressed


  call VTI_readDataset_raw(base64Str,dataType,headerType,compressed, &
                           fileContent,label)
  dataset = as_Int(base64Str,headerType,compressed,dataType)

end function VTI_readDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief Read real dataset from a VTK image data (*.vti) file.
!> @details https://vtk.org/Wiki/VTK_XML_Formats
!--------------------------------------------------------------------------------------------------
function VTI_readDataset_real(fileContent,label) result(dataset)

  character(len=*), intent(in) :: &
    label, &
    fileContent
  real(pREAL),  dimension(:), allocatable :: &
    dataset

  character(len=:), allocatable :: dataType, headerType, base64Str
  logical :: compressed


  call VTI_readDataset_raw(base64Str,dataType,headerType,compressed, &
                           fileContent,label)
  dataset = as_real(base64Str,headerType,compressed,dataType)

end function VTI_readDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Read dataset as raw data (base64 string) from a VTK image data (*.vti) file.
!> @details https://vtk.org/Wiki/VTK_XML_Formats
!--------------------------------------------------------------------------------------------------
subroutine VTI_readDataset_raw(base64Str,dataType,headerType,compressed, &
                               fileContent,label)

  character(len=*), intent(in) :: &
    label, &
    fileContent
  character(len=:), allocatable, intent(out) :: dataType, headerType, base64Str
  logical, intent(out) :: compressed

  logical :: inFile, inImage
  integer(pI64) :: &
    startPos, endPos, &
    s


  inFile = .false.
  inImage = .false.
  startPos = 1_pI64
  outer: do while (startPos < len(fileContent,kind=pI64))
    endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
    if (endPos < startPos) endPos = len(fileContent,kind=pI64)                                      ! end of file without new line

    if (.not. inFile) then
      if (index(fileContent(startPos:endPos),'<VTKFile',kind=pI64) /= 0_pI64) then
        inFile = .true.
        call checkFileFormat(fileContent(startPos:endPos))
        headerType = merge('UInt64','UInt32',getXMLValue(fileContent(startPos:endPos),'header_type')=='UInt64')
        compressed  = getXMLValue(fileContent(startPos:endPos),'compressor') == 'vtkZLibDataCompressor'
      end if
    else
      if (.not. inImage) then
        if (index(fileContent(startPos:endPos),'<ImageData',kind=pI64) /= 0_pI64) then
          inImage = .true.
        end if
      else
        if (index(fileContent(startPos:endPos),'<CellData',kind=pI64) /= 0_pI64) then
          do while (index(fileContent(startPos:endPos),'</CellData>',kind=pI64) == 0_pI64)
            if (index(fileContent(startPos:endPos),'<DataArray',kind=pI64) /= 0_pI64 .and. &
                 getXMLValue(fileContent(startPos:endPos),'Name') == label ) then

              if (getXMLValue(fileContent(startPos:endPos),'format') /= 'binary') &
                call IO_error(844_pI16, 'dataset', label, 'not in binary format', emph = [2])
              dataType = getXMLValue(fileContent(startPos:endPos),'type')

              startPos = endPos + 2_pI64
              endPos  = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
              s = startPos + verify(fileContent(startPos:endPos),IO_WHITESPACE,kind=pI64) -1_pI64   ! start (no leading whitespace)
              base64Str = fileContent(s:endPos)
              exit outer
            end if
            startPos = endPos + 2_pI64
            endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
          end do
        end if
      end if
    end if

    startPos = endPos + 2_pI64

  end do outer

  if (.not. allocated(base64Str)) call IO_error(844_pI16, 'dataset', label, 'not found', emph = [2])

end subroutine VTI_readDataset_raw


!--------------------------------------------------------------------------------------------------
!> @brief Read cells, size, and origin, and cell data labels of an VTK image data (*.vti) file.
!> @details https://vtk.org/Wiki/VTK_XML_Formats
!--------------------------------------------------------------------------------------------------
subroutine VTI_readGeometry(cells,geomSize,origin,labels, &
                            fileContent)

  integer,     dimension(3), intent(out) :: &
    cells                                                                                           ! # of cells (across all processes!)
  real(pREAL), dimension(3), intent(out) :: &
    geomSize, &                                                                                     ! size (across all processes!)
    origin                                                                                          ! origin (across all processes!)
  character(len=pSTRLEN), allocatable, dimension(:), intent(out) :: &
    labels                                                                                          ! cell data labels
  character(len=*), intent(in) :: &
    fileContent

  character(len=:), allocatable :: headerType
  logical :: inFile, inImage, compressed
  integer(pI64) :: &
    startPos, endPos


  cells = -1
  geomSize = -1.0_pREAL

  inFile = .false.
  inImage = .false.
  startPos = 1_pI64

  do while (startPos < len(fileContent,kind=pI64))
    endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
    if (endPos < startPos) endPos = len(fileContent,kind=pI64)                                      ! end of file without new line

    if (.not. inFile) then
      if (index(fileContent(startPos:endPos),'<VTKFile',kind=pI64) /= 0_pI64) then
        inFile = .true.
        call checkFileFormat(fileContent(startPos:endPos))
        headerType = merge('UInt64','UInt32',getXMLValue(fileContent(startPos:endPos),'header_type')=='UInt64')
        compressed = getXMLValue(fileContent(startPos:endPos),'compressor') == 'vtkZLibDataCompressor'
      end if
    else
      if (.not. inImage) then
        if (index(fileContent(startPos:endPos),'<ImageData',kind=pI64) /= 0_pI64) then
          inImage = .true.
          call cellsSizeOrigin(cells,geomSize,origin,fileContent(startPos:endPos))
        end if
      else
        if (index(fileContent(startPos:endPos),'<CellData',kind=pI64) /= 0_pI64) then
          call cell_labels(labels,fileContent(startPos:))
          exit
        end if
      end if
    end if

    startPos = endPos + 2_pI64

  end do

  if (any(geomSize<=0)) call IO_error(844_pI16, 'one or more entries <= 0 for', 'size', emph=[2])
  if (any(cells<1))     call IO_error(844_pI16, 'one or more entries < 1 for', 'cells', emph=[2])

end subroutine VTI_readGeometry


!--------------------------------------------------------------------------------------------------
!> @brief Determine size and origin from coordinates.
!--------------------------------------------------------------------------------------------------
subroutine cellsSizeOrigin(c,s,o,header)

  integer, dimension(3),     intent(out) :: c
  real(pREAL), dimension(3), intent(out) :: s,o
  character(len=*),          intent(in) :: header

  character(len=:), allocatable, dimension(:) :: temp
  real(pREAL), dimension(3) :: delta
  integer :: i


  temp = [getXMLValue(header,'Direction')]
  if (temp(1) /= '1 0 0 0 1 0 0 0 1' .and. temp(1) /= '') &                                         ! https://discourse.vtk.org/t/vti-specification/6526
    call IO_error(844_pI16, 'wrong coordinate order', temp(1), emph=[2])

  call tokenize(getXMLValue(header,'WholeExtent'),' ',temp)
  if (any([(IO_strAsInt(temp(i)),i=1,5,2)] /= 0)) &
    call IO_error(844_pI16, 'coordinate start not at 0', getXMLValue(header,'WholeExtent'), emph=[2])
  c = [(IO_strAsInt(temp(i)),i=2,6,2)]

  call tokenize(getXMLValue(header,'Spacing'),' ',temp)
  delta = [(IO_strAsReal(temp(i)),i=1,3)]
  s = delta * real(c,pREAL)

  call tokenize(getXMLValue(header,'Origin'),' ',temp)
  o = [(IO_strAsReal(temp(i)),i=1,3)]

end subroutine cellsSizeOrigin


!--------------------------------------------------------------------------------------------------
!> @brief Get labels of all cell-based datasets.
!--------------------------------------------------------------------------------------------------
subroutine cell_labels(labels,file_content)

  character(len=pSTRLEN), allocatable, dimension(:), intent(out) :: labels                          !< labels of cell data
  character(len=*), intent(in) :: file_content

  character(len=pSTRLEN) :: label
  integer(pI64) :: startPos, endPos


  startPos = 1_pI64
  endPos = startPos + index(file_content(startPos:),IO_EOL,kind=pI64) - 2_pI64

  allocate(labels(0))

  do while (index(file_content(startPos:endPos),'</CellData>',kind=pI64) == 0_pI64)
    if (index(file_content(startPos:endPos),'<DataArray',kind=pI64) /= 0_pI64) then
      label = getXMLValue(file_content(startPos:endPos),'Name')
      if (any(labels == label)) then
        call IO_error(844_pI16, 'repeated label', trim(label), emph = [2])
      else
        labels = [labels, label]
      end if
    end if
    startPos = endPos + 2_pI64
    endPos = startPos + index(file_content(startPos:),IO_EOL,kind=pI64) - 2_pI64
  end do

end subroutine cell_labels


!--------------------------------------------------------------------------------------------------
!> @brief Interpret Base64 string in vtk XML file as integer of default kind.
!--------------------------------------------------------------------------------------------------
function as_Int(base64Str,headerType,compressed,dataType)

  character(len=*), intent(in) :: base64Str, &                                                      ! base64 encoded string
                                  headerType, &                                                     ! header type (UInt32 or Uint64)
                                  dataType                                                          ! data type (Int32, Int64, Float32, Float64)
  logical,          intent(in) :: compressed                                                        ! indicate whether data is zlib compressed

  integer, dimension(:), allocatable :: as_Int


  select case(dataType)
    case('Int32')
      as_Int = int(prec_bytesToC_INT32_T(asBytes(base64Str,headerType,compressed)))
    case('Int64')
      as_Int = int(prec_bytesToC_INT64_T(asBytes(base64Str,headerType,compressed)))
    case('Float32')
      as_Int = int(prec_bytesToC_FLOAT  (asBytes(base64Str,headerType,compressed)))
    case('Float64')
      as_Int = int(prec_bytesToC_DOUBLE (asBytes(base64Str,headerType,compressed)))
    case default
      call IO_error(844_pI16,'unknown data type',trim(dataType), emph=[2])
  end select

end function as_Int


!--------------------------------------------------------------------------------------------------
!> @brief Interpret Base64 string in vtk XML file as real of kind pREAL.
!--------------------------------------------------------------------------------------------------
function as_real(base64Str,headerType,compressed,dataType)

  character(len=*), intent(in) :: base64Str, &                                                      ! base64 encoded string
                                  headerType, &                                                     ! header type (UInt32 or Uint64)
                                  dataType                                                          ! data type (Int32, Int64, Float32, Float64)
  logical,          intent(in) :: compressed                                                        ! indicate whether data is zlib compressed

  real(pREAL), dimension(:), allocatable :: as_real


  select case(dataType)
    case('Int32')
      as_real = real(prec_bytesToC_INT32_T(asBytes(base64Str,headerType,compressed)),pREAL)
    case('Int64')
      as_real = real(prec_bytesToC_INT64_T(asBytes(base64Str,headerType,compressed)),pREAL)
    case('Float32')
      as_real = real(prec_bytesToC_FLOAT  (asBytes(base64Str,headerType,compressed)),pREAL)
    case('Float64')
      as_real = real(prec_bytesToC_DOUBLE (asBytes(base64Str,headerType,compressed)),pREAL)
    case default
      call IO_error(844_pI16,'unknown data type',trim(dataType), emph=[2])
  end select

end function as_real


!--------------------------------------------------------------------------------------------------
!> @brief Interpret Base64 string in vtk XML file as bytes.
!--------------------------------------------------------------------------------------------------
function asBytes(base64Str,headerType,compressed) result(bytes)

  character(len=*), intent(in) :: base64Str, &                                                      ! base64 encoded string
                                  headerType                                                        ! header type (UInt32 or Uint64)
  logical,          intent(in) :: compressed                                                        ! indicate whether data is zlib compressed

  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes


  if (compressed) then
    bytes = asBytes_compressed(base64Str,headerType)
  else
    bytes = asBytes_uncompressed(base64Str,headerType)
  end if

end function asBytes


!--------------------------------------------------------------------------------------------------
!> @brief Interpret compressed Base64 string in vtk XML file as bytes.
!> @details A compressed Base64 string consists of a header block and a data block
! [#blocks/#u-size/#p-size/#c-size-1/#c-size-2/.../#c-size-#blocks][DATA-1/DATA-2...]
! #blocks = Number of blocks
! #u-size = Block size before compression
! #p-size = Size of last partial block (zero if it not needed)
! #c-size-i = Size in bytes of block i after compression
!--------------------------------------------------------------------------------------------------
function asBytes_compressed(base64Str,headerType) result(bytes_inflated)

  character(len=*), intent(in) :: base64Str, &                                                      ! base64 encoded string
                                  headerType                                                        ! header type (UInt32 or Uint64)
  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes_inflated

  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes_deflated
  integer(pI64), dimension(:), allocatable :: temp, size_inflated, size_deflated
  integer(pI64) :: headerLen, nBlock, b,s,e


  if      (headerType == 'UInt32') then
    temp = int(prec_bytesToC_INT32_T(base64_to_bytes(base64Str(:base64_nChar(4_pI64)))),pI64)
    nBlock = int(temp(1),pI64)
    headerLen = 4_pI64 * (3_pI64 + nBlock)
    temp = int(prec_bytesToC_INT32_T(base64_to_bytes(base64Str(:base64_nChar(headerLen)))),pI64)
  else if (headerType == 'UInt64') then
    temp = int(prec_bytesToC_INT64_T(base64_to_bytes(base64Str(:base64_nChar(8_pI64)))),pI64)
    nBlock = int(temp(1),pI64)
    headerLen = 8_pI64 * (3_pI64 + nBlock)
    temp = int(prec_bytesToC_INT64_T(base64_to_bytes(base64Str(:base64_nChar(headerLen)))),pI64)
  end if

  allocate(size_inflated(nBlock),source=temp(2))
  size_inflated(nBlock) = merge(temp(3),temp(2),temp(3)/=0_pI64)
  size_deflated = temp(4:)
  bytes_deflated = base64_to_bytes(base64Str(base64_nChar(headerLen)+1_pI64:))

  allocate(bytes_inflated(sum(size_inflated)))
  e = 0_pI64
  do b = 1, nBlock
    s = e + 1_pI64
    e = s + size_deflated(b) - 1_pI64
    bytes_inflated(sum(size_inflated(:b-1))+1_pI64:sum(size_inflated(:b))) = zlib_inflate(bytes_deflated(s:e),size_inflated(b))
  end do

end function asBytes_compressed


!--------------------------------------------------------------------------------------------------
!> @brief Interprete uncompressed Base64 string in vtk XML file as bytes.
!> @details An uncompressed Base64 string consists of N headers blocks and a N data blocks
![#bytes-1/DATA-1][#bytes-2/DATA-2]...
!--------------------------------------------------------------------------------------------------
function asBytes_uncompressed(base64Str,headerType) result(bytes)

  character(len=*), intent(in) :: base64Str, &                                                      ! base64 encoded string
                                  headerType                                                        ! header type (UInt32 or Uint64)
  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes

  integer(pI64) :: s
  integer(pI64), dimension(1) :: nByte


  allocate(bytes(0))

  s=0_pI64
  if      (headerType == 'UInt32') then
    do while(s+base64_nChar(4_pI64)<(len(base64Str,pI64)))
      nByte = int(prec_bytesToC_INT32_T(base64_to_bytes(base64Str(s+1_pI64:s+base64_nChar(4_pI64)))),pI64)
      bytes = [bytes,base64_to_bytes(base64Str(s+1_pI64:s+base64_nChar(4_pI64+nByte(1))),5_pI64)]
      s = s + base64_nChar(4_pI64+nByte(1))
    end do
  else if (headerType == 'UInt64') then
    do while(s+base64_nChar(8_pI64)<(len(base64Str,pI64)))
      nByte = int(prec_bytesToC_INT64_T(base64_to_bytes(base64Str(s+1_pI64:s+base64_nChar(8_pI64)))),pI64)
      bytes = [bytes,base64_to_bytes(base64Str(s+1_pI64:s+base64_nChar(8_pI64+nByte(1))),9_pI64)]
      s = s + base64_nChar(8_pI64+nByte(1))
    end do
  end if

end function asBytes_uncompressed


!--------------------------------------------------------------------------------------------------
!> @brief Get XML string value for given key.
!--------------------------------------------------------------------------------------------------
pure function getXMLValue(line,key)

  character(len=*), intent(in)  :: line, key
  character(len=:), allocatable :: getXMLValue

  integer :: s,e
#ifdef __INTEL_COMPILER
  character :: q
#endif


  s = index(line," "//key,back=.true.)
  if (s==0) then
    getXMLValue = ''
  else
    e = s + 1 + scan(line(s+1:),"'"//'"')
    if (scan(line(s:e-2),'=') == 0) then
      getXMLValue = ''
    else
      s = e
!https://community.intel.com/t5/Intel-Fortran-Compiler/ICE-for-merge-with-strings/m-p/1207204#M151657
#ifdef __INTEL_COMPILER
      q = line(s-1:s-1)
      e = s + index(line(s:),q) - 1
#else
      e = s + index(line(s:),merge("'",'"',line(s-1:s-1)=="'")) - 1
#endif
      getXMLValue = line(s:e-1)
    end if
  end if

end function getXMLValue


!--------------------------------------------------------------------------------------------------
!> @brief Check for supported file format variants.
!--------------------------------------------------------------------------------------------------
subroutine checkFileFormat(line)

  character(len=*),intent(in) :: line

  character(len=:), allocatable :: val


  val = getXMLValue(line,'type')
  if (val /= 'ImageData') &
    call IO_error(844_pI16, 'type', val, 'is not', 'ImageData',emph=[2,4])

  val = getXMLValue(line,'byte_order')
  if (val /= 'LittleEndian') &
    call IO_error(844_pI16, 'byte_order', val, 'is not', 'LittleEndian',emph=[2,4])

  val = getXMLValue(line,'compressor')
  if (val /= '' .and. val /= 'vtkZLibDataCompressor') &
    call IO_error(844_pI16, 'compressor', val, 'is not', 'vtkZLibDataCompressor',emph=[2,4])

end subroutine checkFileFormat

end module VTI
