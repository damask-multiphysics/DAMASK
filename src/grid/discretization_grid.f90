!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse geometry file to set up discretization and geometry for nonlocal model
!--------------------------------------------------------------------------------------------------
module discretization_grid
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  use prec
  use parallelization
  use system_routines
  use base64
  use zlib
  use DAMASK_interface
  use IO
  use config
  use results
  use discretization
  use geometry_plastic_nonlocal

  implicit none
  private

  integer,     dimension(3), public, protected :: &
    cells                                                                                           !< (global) cells
  integer,                   public, protected :: &
    cells3, &                                                                                       !< (local) cells in 3rd direction
    cells3Offset                                                                                    !< (local) cells offset in 3rd direction
  real(pReal), dimension(3), public, protected :: &
    geomSize                                                                                        !< (global) physical size
  real(pReal),               public, protected :: &
    size3, &                                                                                        !< (local) size in 3rd direction
    size3offset                                                                                     !< (local) size offset in 3rd direction

  public :: &
    discretization_grid_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads the geometry file to obtain information on discretization
!--------------------------------------------------------------------------------------------------
subroutine discretization_grid_init(restart)

  logical, intent(in) :: restart

  include 'fftw3-mpi.f03'
  real(pReal), dimension(3) :: &
    mySize, &                                                                                       !< domain size of this process
    origin                                                                                          !< (global) distance to origin
  integer,     dimension(3) :: &
    myGrid                                                                                          !< domain cells of this process

  integer,     dimension(:),   allocatable :: &
    materialAt, materialAt_global

  integer :: &
    j, &
    debug_element, debug_ip
  integer(MPI_INTEGER_KIND) :: err_MPI
  integer(C_INTPTR_T) :: &
    devNull, z, z_offset
  integer, dimension(worldsize) :: &
    displs, sendcounts
  character(len=:), allocatable :: &
    fileContent, fname


  print'(/,1x,a)', '<<<+-  discretization_grid init  -+>>>'; flush(IO_STDOUT)


  if (worldrank == 0) then
    fileContent = IO_read(interface_geomFile)
    call readVTI(cells,geomSize,origin,materialAt_global,fileContent)
    fname = interface_geomFile
    if (scan(fname,'/') /= 0) fname = fname(scan(fname,'/',.true.)+1:)
    call results_openJobFile(parallel=.false.)
    call results_writeDataset_str(fileContent,'setup',fname,'geometry definition (grid solver)')
    call results_closeJobFile
  else
    allocate(materialAt_global(0))                                                                  ! needed for IntelMPI
  end if


  call MPI_Bcast(cells,3_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD, err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  if (cells(1) < 2) call IO_error(844, ext_msg='cells(1) must be larger than 1')
  call MPI_Bcast(geomSize,3_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD, err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call MPI_Bcast(origin,3_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD, err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  print'(/,1x,a,i0,a,i0,a,i0)',            'cells:  ', cells(1),    ' × ', cells(2),    ' × ', cells(3)
  print  '(1x,a,es8.2,a,es8.2,a,es8.2,a)', 'size:   ', geomSize(1), ' × ', geomSize(2), ' × ', geomSize(3), ' m³'
  print  '(1x,a,es8.2,a,es8.2,a,es8.2,a)', 'origin: ', origin(1),   ' ',   origin(2),   ' ',   origin(3), ' m'

  if (worldsize>cells(3)) call IO_error(894, ext_msg='number of processes exceeds cells(3)')

  call fftw_mpi_init
  devNull = fftw_mpi_local_size_3d(int(cells(3),C_INTPTR_T), &
                                   int(cells(2),C_INTPTR_T), &
                                   int(cells(1),C_INTPTR_T)/2+1, &
                                   PETSC_COMM_WORLD, &
                                   z, &                                                             ! domain cells size along z
                                   z_offset)                                                        ! domain cells offset along z
  if (z==0_C_INTPTR_T) call IO_error(894, ext_msg='Cannot distribute MPI processes')

  cells3       = int(z)
  cells3Offset = int(z_offset)
  size3       = geomSize(3)*real(cells3,pReal)      /real(cells(3),pReal)
  size3Offset = geomSize(3)*real(cells3Offset,pReal)/real(cells(3),pReal)
  myGrid = [cells(1:2),cells3]
  mySize = [geomSize(1:2),size3]

  call MPI_Gather(product(cells(1:2))*cells3Offset, 1_MPI_INTEGER_KIND,MPI_INTEGER,displs,&
                  1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call MPI_Gather(product(myGrid),               1_MPI_INTEGER_KIND,MPI_INTEGER,sendcounts,&
                  1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  allocate(materialAt(product(myGrid)))
  call MPI_Scatterv(materialAt_global,sendcounts,displs,MPI_INTEGER,materialAt,size(materialAt),&
                    MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  call discretization_init(materialAt, &
                           IPcoordinates0(myGrid,mySize,cells3Offset), &
                           Nodes0(myGrid,mySize,cells3Offset),&
                           merge((cells(1)+1) * (cells(2)+1) * (cells3+1),&                         ! write top layer...
                                 (cells(1)+1) * (cells(2)+1) *  cells3,&                            ! ...unless not last process
                                 worldrank+1==worldsize))

!--------------------------------------------------------------------------------------------------
! store geometry information for post processing
  if (.not. restart) then
    call results_openJobFile
    call results_closeGroup(results_addGroup('geometry'))
    call results_addAttribute('cells', cells,   '/geometry')
    call results_addAttribute('size',  geomSize,'/geometry')
    call results_addAttribute('origin',origin,  '/geometry')
    call results_closeJobFile
  end if

!--------------------------------------------------------------------------------------------------
! geometry information required by the nonlocal CP model
  call geometry_plastic_nonlocal_setIPvolume(reshape([(product(mySize/real(myGrid,pReal)),j=1,product(myGrid))], &
                                                     [1,product(myGrid)]))
  call geometry_plastic_nonlocal_setIParea        (cellSurfaceArea(mySize,myGrid))
  call geometry_plastic_nonlocal_setIPareaNormal  (cellSurfaceNormal(product(myGrid)))
  call geometry_plastic_nonlocal_setIPneighborhood(IPneighborhood(myGrid))

!-------------------------------------------------------------------------------------------------
! debug parameters
  debug_element = config_debug%get_asInt('element',defaultVal=1)
  if (debug_element < 1 .or. debug_element > product(myGrid)) call IO_error(602,ext_msg='element')
  debug_ip      = config_debug%get_asInt('integrationpoint',defaultVal=1)
  if (debug_ip /= 1)                                          call IO_error(602,ext_msg='IP')

end subroutine discretization_grid_init


!--------------------------------------------------------------------------------------------------
!> @brief Parse vtk image data (.vti)
!> @details https://vtk.org/Wiki/VTK_XML_Formats
!--------------------------------------------------------------------------------------------------
subroutine readVTI(cells,geomSize,origin,material, &
                   fileContent)

  integer,     dimension(3), intent(out) :: &
    cells                                                                                           ! cells (across all processes!)
  real(pReal), dimension(3), intent(out) :: &
    geomSize, &                                                                                     ! size (across all processes!)
    origin                                                                                          ! origin (across all processes!)
  integer,     dimension(:), intent(out), allocatable :: &
    material
  character(len=*),          intent(in) :: &
    fileContent

  character(len=:), allocatable :: dataType, headerType
  logical :: inFile,inImage,gotCellData,compressed
  integer(pI64) :: &
    startPos, endPos, &
    s


  cells = -1
  geomSize = -1.0_pReal

  inFile         = .false.
  inImage        = .false.
  gotCelldata    = .false.

!--------------------------------------------------------------------------------------------------
! parse XML file
  startPos = 1_pI64
  do while (startPos < len(fileContent,kind=pI64))
    endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
    if (endPos < startPos) endPos = len(fileContent,kind=pI64)                                      ! end of file without new line

    if (.not. inFile) then
      if (index(fileContent(startPos:endPos),'<VTKFile',kind=pI64) /= 0_pI64) then
        inFile = .true.
        if (.not. fileFormatOk(fileContent(startPos:endPos))) call IO_error(error_ID = 844, ext_msg='file format')
        headerType = merge('UInt64','UInt32',getXMLValue(fileContent(startPos:endPos),'header_type')=='UInt64')
        compressed  = getXMLValue(fileContent(startPos:endPos),'compressor') == 'vtkZLibDataCompressor'
      end if
    else
      if (.not. inImage) then
        if (index(fileContent(startPos:endPos),'<ImageData',kind=pI64) /= 0_pI64) then
          inImage = .true.
          call cellsSizeOrigin(cells,geomSize,origin,fileContent(startPos:endPos))
        end if
      else
        if (index(fileContent(startPos:endPos),'<CellData',kind=pI64) /= 0_pI64) then
          gotCellData = .true.
          do while (index(fileContent(startPos:endPos),'</CellData>',kind=pI64) == 0_pI64)
            if (index(fileContent(startPos:endPos),'<DataArray',kind=pI64) /= 0_pI64 .and. &
                 getXMLValue(fileContent(startPos:endPos),'Name') == 'material' ) then

              if (getXMLValue(fileContent(startPos:endPos),'format') /= 'binary') &
                call IO_error(error_ID = 844, ext_msg='format (material)')
              dataType = getXMLValue(fileContent(startPos:endPos),'type')

              startPos = endPos + 2_pI64
              endPos  = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
              s = startPos + verify(fileContent(startPos:endPos),IO_WHITESPACE,kind=pI64) -1_pI64   ! start (no leading whitespace)
              material = as_Int(fileContent(s:endPos),headerType,compressed,dataType)
              exit
            end if
            startPos = endPos + 2_pI64
            endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
          end do
        end if
      end if
    end if

    if (gotCellData) exit
    startPos = endPos + 2_pI64

  end do

  if (.not. allocated(material))        call IO_error(error_ID = 844, ext_msg='material data not found')
  if (size(material) /= product(cells)) call IO_error(error_ID = 844, ext_msg='size(material)')
  if (any(geomSize<=0))                 call IO_error(error_ID = 844, ext_msg='size')
  if (any(cells<1))                     call IO_error(error_ID = 844, ext_msg='cells')
  material = material + 1
  if (any(material<1))                 call IO_error(error_ID = 844, ext_msg='material ID < 0')

  contains

  !------------------------------------------------------------------------------------------------
  !> @brief determine size and origin from coordinates
  !------------------------------------------------------------------------------------------------
  subroutine cellsSizeOrigin(c,s,o,header)

    integer, dimension(3),     intent(out) :: c
    real(pReal), dimension(3), intent(out) :: s,o
    character(len=*),          intent(in) :: header

    character(len=:), allocatable :: temp
    real(pReal), dimension(:), allocatable :: delta
    integer :: i


    temp = getXMLValue(header,'Direction')
    if (temp /= '1 0 0 0 1 0 0 0 1' .and. temp /= '') &                                             ! https://discourse.vtk.org/t/vti-specification/6526
      call IO_error(error_ID = 844, ext_msg = 'coordinate order')

    temp = getXMLValue(header,'WholeExtent')
    if (any([(IO_intValue(temp,IO_stringPos(temp),i),i=1,5,2)] /= 0)) &
      call IO_error(error_ID = 844, ext_msg = 'coordinate start')
    c = [(IO_intValue(temp,IO_stringPos(temp),i),i=2,6,2)]

    temp = getXMLValue(header,'Spacing')
    delta = [(IO_floatValue(temp,IO_stringPos(temp),i),i=1,3)]
    s = delta * real(c,pReal)

    temp = getXMLValue(header,'Origin')
    o = [(IO_floatValue(temp,IO_stringPos(temp),i),i=1,3)]

  end subroutine


  !------------------------------------------------------------------------------------------------
  !> @brief Interpret Base64 string in vtk XML file as integer of default kind
  !------------------------------------------------------------------------------------------------
  function as_Int(base64_str,headerType,compressed,dataType)

    character(len=*), intent(in) :: base64_str, &                                                   ! base64 encoded string
                                    headerType, &                                                   ! header type (UInt32 or Uint64)
                                    dataType                                                        ! data type (Int32, Int64, Float32, Float64)
    logical,          intent(in) :: compressed                                                      ! indicate whether data is zlib compressed

    integer, dimension(:), allocatable :: as_Int

    select case(dataType)
      case('Int32')
        as_Int = int(prec_bytesToC_INT32_T(asBytes(base64_str,headerType,compressed)))
      case('Int64')
        as_Int = int(prec_bytesToC_INT64_T(asBytes(base64_str,headerType,compressed)))
      case('Float32')
        as_Int = int(prec_bytesToC_FLOAT  (asBytes(base64_str,headerType,compressed)))
      case('Float64')
        as_Int = int(prec_bytesToC_DOUBLE (asBytes(base64_str,headerType,compressed)))
      case default
        call IO_error(844,ext_msg='unknown data type: '//trim(dataType))
    end select

  end function as_Int


  !------------------------------------------------------------------------------------------------
  !> @brief Interpret Base64 string in vtk XML file as integer of pReal kind
  !------------------------------------------------------------------------------------------------
  function as_pReal(base64_str,headerType,compressed,dataType)

    character(len=*), intent(in) :: base64_str, &                                                   ! base64 encoded string
                                    headerType, &                                                   ! header type (UInt32 or Uint64)
                                    dataType                                                        ! data type (Int32, Int64, Float32, Float64)
    logical,          intent(in) :: compressed                                                      ! indicate whether data is zlib compressed

    real(pReal), dimension(:), allocatable :: as_pReal

    select case(dataType)
      case('Int32')
        as_pReal = real(prec_bytesToC_INT32_T(asBytes(base64_str,headerType,compressed)),pReal)
      case('Int64')
        as_pReal = real(prec_bytesToC_INT64_T(asBytes(base64_str,headerType,compressed)),pReal)
      case('Float32')
        as_pReal = real(prec_bytesToC_FLOAT  (asBytes(base64_str,headerType,compressed)),pReal)
      case('Float64')
        as_pReal = real(prec_bytesToC_DOUBLE (asBytes(base64_str,headerType,compressed)),pReal)
      case default
        call IO_error(844,ext_msg='unknown data type: '//trim(dataType))
    end select

  end function as_pReal


  !------------------------------------------------------------------------------------------------
  !> @brief Interpret Base64 string in vtk XML file as bytes
  !------------------------------------------------------------------------------------------------
  function asBytes(base64_str,headerType,compressed) result(bytes)

    character(len=*), intent(in) :: base64_str, &                                                   ! base64 encoded string
                                    headerType                                                      ! header type (UInt32 or Uint64)
    logical,          intent(in) :: compressed                                                      ! indicate whether data is zlib compressed

    integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes

    if (compressed) then
      bytes = asBytes_compressed(base64_str,headerType)
    else
      bytes = asBytes_uncompressed(base64_str,headerType)
    end if

  end function asBytes

  !------------------------------------------------------------------------------------------------
  !> @brief Interpret compressed Base64 string in vtk XML file as bytes
  !> @details A compressed Base64 string consists of a header block and a data block
  ! [#blocks/#u-size/#p-size/#c-size-1/#c-size-2/.../#c-size-#blocks][DATA-1/DATA-2...]
  ! #blocks = Number of blocks
  ! #u-size = Block size before compression
  ! #p-size = Size of last partial block (zero if it not needed)
  ! #c-size-i = Size in bytes of block i after compression
  !------------------------------------------------------------------------------------------------
  function asBytes_compressed(base64_str,headerType) result(bytes)

    character(len=*), intent(in) :: base64_str, &                                                   ! base64 encoded string
                                    headerType                                                      ! header type (UInt32 or Uint64)

    integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes, bytes_inflated

    integer(pI64), dimension(:), allocatable :: temp, size_inflated, size_deflated
    integer(pI64) :: headerLen, nBlock, b,s,e


    if      (headerType == 'UInt32') then
      temp = int(prec_bytesToC_INT32_T(base64_to_bytes(base64_str(:base64_nChar(4_pI64)))),pI64)
      nBlock = int(temp(1),pI64)
      headerLen = 4_pI64 * (3_pI64 + nBlock)
      temp = int(prec_bytesToC_INT32_T(base64_to_bytes(base64_str(:base64_nChar(headerLen)))),pI64)
    else if (headerType == 'UInt64') then
      temp = int(prec_bytesToC_INT64_T(base64_to_bytes(base64_str(:base64_nChar(8_pI64)))),pI64)
      nBlock = int(temp(1),pI64)
      headerLen = 8_pI64 * (3_pI64 + nBlock)
      temp = int(prec_bytesToC_INT64_T(base64_to_bytes(base64_str(:base64_nChar(headerLen)))),pI64)
    end if

    allocate(size_inflated(nBlock),source=temp(2))
    size_inflated(nBlock) = merge(temp(3),temp(2),temp(3)/=0_pI64)
    size_deflated = temp(4:)
    bytes_inflated = base64_to_bytes(base64_str(base64_nChar(headerLen)+1_pI64:))

    allocate(bytes(sum(size_inflated)))
    e = 0_pI64
    do b = 1, nBlock
      s = e + 1_pI64
      e = s + size_deflated(b) - 1_pI64
      bytes(sum(size_inflated(:b-1))+1_pI64:sum(size_inflated(:b))) = zlib_inflate(bytes_inflated(s:e),size_inflated(b))
    end do

  end function asBytes_compressed


  !------------------------------------------------------------------------------------------------
  !> @brief Interprete uncompressed Base64 string in vtk XML file as bytes
  !> @details An uncompressed Base64 string consists of N headers blocks and a N data blocks
  ![#bytes-1/DATA-1][#bytes-2/DATA-2]...
  !------------------------------------------------------------------------------------------------
  function asBytes_uncompressed(base64_str,headerType) result(bytes)

    character(len=*), intent(in) :: base64_str, &                                                   ! base64 encoded string
                                    headerType                                                      ! header type (UInt32 or Uint64)

    integer(pI64) :: s
    integer(pI64), dimension(1) :: nByte

    integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes
    allocate(bytes(0))

    s=0_pI64
    if      (headerType == 'UInt32') then
      do while(s+base64_nChar(4_pI64)<(len(base64_str,pI64)))
        nByte = int(prec_bytesToC_INT32_T(base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(4_pI64)))),pI64)
        bytes = [bytes,base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(4_pI64+nByte(1))),5_pI64)]
        s = s + base64_nChar(4_pI64+nByte(1))
      end do
    else if (headerType == 'UInt64') then
      do while(s+base64_nChar(8_pI64)<(len(base64_str,pI64)))
        nByte = int(prec_bytesToC_INT64_T(base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(8_pI64)))),pI64)
        bytes = [bytes,base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(8_pI64+nByte(1))),9_pI64)]
        s = s + base64_nChar(8_pI64+nByte(1))
      end do
    end if

  end function asBytes_uncompressed

  !------------------------------------------------------------------------------------------------
  !> @brief Get XML string value for given key
  !------------------------------------------------------------------------------------------------
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
! https://community.intel.com/t5/Intel-Fortran-Compiler/ICE-for-merge-with-strings/m-p/1207204#M151657
#ifdef __INTEL_COMPILER
        q = line(s-1:s-1)
        e = s + index(line(s:),q) - 1
#else
        e = s + index(line(s:),merge("'",'"',line(s-1:s-1)=="'")) - 1
#endif
        getXMLValue = line(s:e-1)
      end if
    end if

  end function


  !------------------------------------------------------------------------------------------------
  !> @brief check for supported file format
  !------------------------------------------------------------------------------------------------
  pure function fileFormatOk(line)

    character(len=*),intent(in) :: line
    logical :: fileFormatOk

    fileFormatOk = getXMLValue(line,'type')       == 'ImageData' .and. &
                   getXMLValue(line,'byte_order') == 'LittleEndian' .and. &
                   getXMLValue(line,'compressor') /= 'vtkLZ4DataCompressor' .and. &
                   getXMLValue(line,'compressor') /= 'vtkLZMADataCompressor'

  end function fileFormatOk

end subroutine readVTI


!---------------------------------------------------------------------------------------------------
!> @brief Calculate undeformed position of IPs/cell centers (pretend to be an element)
!---------------------------------------------------------------------------------------------------
function IPcoordinates0(cells,geomSize,cells3Offset)

  integer,     dimension(3), intent(in) :: cells                                                    ! cells (for this process!)
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,                   intent(in) :: cells3Offset                                             ! cells(3) offset

  real(pReal), dimension(3,product(cells))  :: ipCoordinates0

  integer :: &
    a,b,c, &
    i


  i = 0
  do c = 1, cells(3); do b = 1, cells(2); do a = 1, cells(1)
    i = i + 1
    IPcoordinates0(1:3,i) = geomSize/real(cells,pReal) * (real([a,b,cells3Offset+c],pReal) -0.5_pReal)
  end do; end do; end do

end function IPcoordinates0


!---------------------------------------------------------------------------------------------------
!> @brief Calculate position of undeformed nodes (pretend to be an element)
!---------------------------------------------------------------------------------------------------
pure function nodes0(cells,geomSize,cells3Offset)

  integer,     dimension(3), intent(in) :: cells                                                    ! cells (for this process!)
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,                   intent(in) :: cells3Offset                                             ! cells(3) offset

  real(pReal), dimension(3,product(cells+1)) :: nodes0

  integer :: &
    a,b,c, &
    n

  n = 0
  do c = 0, cells3; do b = 0, cells(2); do a = 0, cells(1)
    n = n + 1
    nodes0(1:3,n) = geomSize/real(cells,pReal) * real([a,b,cells3Offset+c],pReal)
  end do; end do; end do

end function nodes0


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP interface areas
!--------------------------------------------------------------------------------------------------
pure function cellSurfaceArea(geomSize,cells)

  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,     dimension(3), intent(in) :: cells                                                    ! cells (for this process!)

  real(pReal), dimension(6,1,product(cells)) :: cellSurfaceArea


  cellSurfaceArea(1:2,1,:) = geomSize(2)/real(cells(2)) * geomSize(3)/real(cells(3))
  cellSurfaceArea(3:4,1,:) = geomSize(3)/real(cells(3)) * geomSize(1)/real(cells(1))
  cellSurfaceArea(5:6,1,:) = geomSize(1)/real(cells(1)) * geomSize(2)/real(cells(2))

end function cellSurfaceArea


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP interface areas normals
!--------------------------------------------------------------------------------------------------
pure function cellSurfaceNormal(nElems)

  integer, intent(in) :: nElems

  real(pReal), dimension(3,6,1,nElems) :: cellSurfaceNormal

  cellSurfaceNormal(1:3,1,1,:) = spread([+1.0_pReal, 0.0_pReal, 0.0_pReal],2,nElems)
  cellSurfaceNormal(1:3,2,1,:) = spread([-1.0_pReal, 0.0_pReal, 0.0_pReal],2,nElems)
  cellSurfaceNormal(1:3,3,1,:) = spread([ 0.0_pReal,+1.0_pReal, 0.0_pReal],2,nElems)
  cellSurfaceNormal(1:3,4,1,:) = spread([ 0.0_pReal,-1.0_pReal, 0.0_pReal],2,nElems)
  cellSurfaceNormal(1:3,5,1,:) = spread([ 0.0_pReal, 0.0_pReal,+1.0_pReal],2,nElems)
  cellSurfaceNormal(1:3,6,1,:) = spread([ 0.0_pReal, 0.0_pReal,-1.0_pReal],2,nElems)

end function cellSurfaceNormal


!--------------------------------------------------------------------------------------------------
!> @brief Build IP neighborhood relations
!--------------------------------------------------------------------------------------------------
pure function IPneighborhood(cells)

  integer, dimension(3), intent(in) :: cells                                                        ! cells (for this process!)

  integer, dimension(3,6,1,product(cells)) :: IPneighborhood                                        !< 6 neighboring IPs as [element ID, IP ID, face ID]

  integer :: &
   x,y,z, &
   e

  e = 0
  do z = 0,cells(3)-1; do y = 0,cells(2)-1; do x = 0,cells(1)-1
    e = e + 1
    ! element ID
    IPneighborhood(1,1,1,e) = z * cells(1) * cells(2) &
                            + y * cells(1) &
                            + modulo(x+1,cells(1)) &
                            + 1
    IPneighborhood(1,2,1,e) = z * cells(1) * cells(2) &
                            + y * cells(1) &
                            + modulo(x-1,cells(1)) &
                            + 1
    IPneighborhood(1,3,1,e) = z * cells(1) * cells(2) &
                            + modulo(y+1,cells(2)) * cells(1) &
                            + x &
                            + 1
    IPneighborhood(1,4,1,e) = z * cells(1) * cells(2) &
                            + modulo(y-1,cells(2)) * cells(1) &
                            + x &
                            + 1
    IPneighborhood(1,5,1,e) = modulo(z+1,cells(3)) * cells(1) * cells(2) &
                            + y * cells(1) &
                            + x &
                            + 1
    IPneighborhood(1,6,1,e) = modulo(z-1,cells(3)) * cells(1) * cells(2) &
                            + y * cells(1) &
                            + x &
                            + 1
    ! IP ID
    IPneighborhood(2,:,1,e) = 1

    ! face ID
    IPneighborhood(3,1,1,e) = 2
    IPneighborhood(3,2,1,e) = 1
    IPneighborhood(3,3,1,e) = 4
    IPneighborhood(3,4,1,e) = 3
    IPneighborhood(3,5,1,e) = 6
    IPneighborhood(3,6,1,e) = 5

  end do; end do; end do

end function IPneighborhood


end module discretization_grid
