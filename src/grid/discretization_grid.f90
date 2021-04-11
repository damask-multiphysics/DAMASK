!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse geometry file to set up discretization and geometry for nonlocal model
!--------------------------------------------------------------------------------------------------
module discretization_grid
#include <petsc/finclude/petscsys.h>
  use PETScsys

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
    grid                                                                                            !< (global) grid
  integer,                   public, protected :: &
    grid3, &                                                                                        !< (local) grid in 3rd direction
    grid3Offset                                                                                     !< (local) grid offset in 3rd direction
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
    myGrid                                                                                          !< domain grid of this process

  integer,     dimension(:),   allocatable :: &
    materialAt, materialAt_global

  integer :: &
    j, &
    debug_element, debug_ip, &
    ierr
  integer(C_INTPTR_T) :: &
    devNull, z, z_offset
  integer, dimension(worldsize) :: &
    displs, sendcounts

  print'(/,a)', ' <<<+-  discretization_grid init  -+>>>'; flush(IO_STDOUT)

  if(worldrank == 0) then
    call readVTR(grid,geomSize,origin,materialAt_global)
  else
    allocate(materialAt_global(0))                                                                  ! needed for IntelMPI
  endif


  call MPI_Bcast(grid,3,MPI_INTEGER,0,PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'
  if (grid(1) < 2) call IO_error(844, ext_msg='cells(1) must be larger than 1')
  call MPI_Bcast(geomSize,3,MPI_DOUBLE,0,PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'
  call MPI_Bcast(origin,3,MPI_DOUBLE,0,PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'

  print'(/,a,3(i12  ))',  ' cells    a b c: ', grid
  print'(a,3(es12.5))',   ' size     x y z: ', geomSize
  print'(a,3(es12.5))',   ' origin   x y z: ', origin

  if(worldsize>grid(3)) call IO_error(894, ext_msg='number of processes exceeds grid(3)')

  call fftw_mpi_init
  devNull = fftw_mpi_local_size_3d(int(grid(3),C_INTPTR_T), &
                                   int(grid(2),C_INTPTR_T), &
                                   int(grid(1),C_INTPTR_T)/2+1, &
                                   PETSC_COMM_WORLD, &
                                   z, &                                                             ! domain grid size along z
                                   z_offset)                                                        ! domain grid offset along z
  if(z==0_C_INTPTR_T) call IO_error(894, ext_msg='Cannot distribute MPI processes')

  grid3       = int(z)
  grid3Offset = int(z_offset)
  size3       = geomSize(3)*real(grid3,pReal)      /real(grid(3),pReal)
  size3Offset = geomSize(3)*real(grid3Offset,pReal)/real(grid(3),pReal)
  myGrid = [grid(1:2),grid3]
  mySize = [geomSize(1:2),size3]

  call MPI_Gather(product(grid(1:2))*grid3Offset,1,MPI_INTEGER,displs,    1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  if (ierr /= 0) error stop 'MPI error'
  call MPI_Gather(product(myGrid),               1,MPI_INTEGER,sendcounts,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  if (ierr /= 0) error stop 'MPI error'

  allocate(materialAt(product(myGrid)))
  call MPI_scatterv(materialAt_global,sendcounts,displs,MPI_INTEGER,materialAt,size(materialAt),MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  if (ierr /= 0) error stop 'MPI error'

  call discretization_init(materialAt, &
                           IPcoordinates0(myGrid,mySize,grid3Offset), &
                           Nodes0(myGrid,mySize,grid3Offset),&
                           merge((grid(1)+1) * (grid(2)+1) * (grid3+1),&                            ! write top layer...
                                 (grid(1)+1) * (grid(2)+1) *  grid3,&                               ! ...unless not last process
                                 worldrank+1==worldsize))

!--------------------------------------------------------------------------------------------------
! store geometry information for post processing
  if(.not. restart) then
    call results_openJobFile
    call results_closeGroup(results_addGroup('geometry'))
    call results_addAttribute('cells', grid,    '/geometry')
    call results_addAttribute('size',  geomSize,'/geometry')
    call results_addAttribute('origin',origin,  '/geometry')
    call results_closeJobFile
  endif

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
!> @brief Parse vtk rectilinear grid (.vtr)
!> @details https://vtk.org/Wiki/VTK_XML_Formats
!--------------------------------------------------------------------------------------------------
subroutine readVTR(grid,geomSize,origin,material)

  integer,     dimension(3), intent(out) :: &
    grid                                                                                            ! grid   (across all processes!)
  real(pReal), dimension(3), intent(out) :: &
    geomSize, &                                                                                     ! size   (across all processes!)
    origin                                                                                          ! origin (across all processes!)
  integer,     dimension(:), intent(out), allocatable :: &
    material

  character(len=:), allocatable :: fileContent, dataType, headerType
  logical :: inFile,inGrid,gotCoordinates,gotCellData,compressed
  integer :: fileUnit, myStat, coord
  integer(pI64) :: &
    fileLength, &                                                                                   !< length of the geom file (in characters)
    startPos, endPos, &
    s

  grid = -1
  geomSize = -1.0_pReal

!--------------------------------------------------------------------------------------------------
! read raw data as stream
  inquire(file = trim(interface_geomFile), size=fileLength)
  open(newunit=fileUnit, file=trim(interface_geomFile), access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if(myStat /= 0) call IO_error(100,ext_msg=trim(interface_geomFile))
  allocate(character(len=fileLength)::fileContent)
  read(fileUnit) fileContent
  close(fileUnit)

  inFile         = .false.
  inGrid         = .false.
  gotCoordinates = .false.
  gotCelldata    = .false.

!--------------------------------------------------------------------------------------------------
! interpret XML file
  startPos = 1_pI64
  do while (startPos < len(fileContent,kind=pI64))
    endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
    if (endPos < startPos) endPos = len(fileContent,kind=pI64)                                      ! end of file without new line

    if(.not. inFile) then
      if(index(fileContent(startPos:endPos),'<VTKFile',kind=pI64) /= 0_pI64) then
        inFile = .true.
        if(.not. fileFormatOk(fileContent(startPos:endPos))) call IO_error(error_ID = 844, ext_msg='file format')
        headerType = merge('UInt64','UInt32',getXMLValue(fileContent(startPos:endPos),'header_type')=='UInt64')
        compressed  = getXMLValue(fileContent(startPos:endPos),'compressor') == 'vtkZLibDataCompressor'
      endif
    else
      if(.not. inGrid) then
        if(index(fileContent(startPos:endPos),'<RectilinearGrid',kind=pI64) /= 0_pI64) inGrid = .true.
      else
        if(index(fileContent(startPos:endPos),'<CellData>',kind=pI64) /= 0_pI64) then
          gotCellData = .true.
          do while (index(fileContent(startPos:endPos),'</CellData>',kind=pI64) == 0_pI64)
            if(index(fileContent(startPos:endPos),'<DataArray',kind=pI64) /= 0_pI64 .and. &
                 getXMLValue(fileContent(startPos:endPos),'Name') == 'material' ) then

              if(getXMLValue(fileContent(startPos:endPos),'format') /= 'binary') &
                call IO_error(error_ID = 844, ext_msg='format (materialpoint)')
              dataType = getXMLValue(fileContent(startPos:endPos),'type')

              startPos = endPos + 2_pI64
              endPos  = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
              s = startPos + verify(fileContent(startPos:endPos),IO_WHITESPACE,kind=pI64) -1_pI64   ! start (no leading whitespace)
              material = as_Int(fileContent(s:endPos),headerType,compressed,dataType)
              exit
            endif
            startPos = endPos + 2_pI64
            endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
          enddo
        elseif(index(fileContent(startPos:endPos),'<Coordinates>',kind=pI64) /= 0_pI64) then
          gotCoordinates = .true.
          startPos = endPos + 2_pI64

          coord = 0
          do while (startPos<fileLength)
            endPos = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
            if(index(fileContent(startPos:endPos),'<DataArray',kind=pI64) /= 0_pI64) then

              if(getXMLValue(fileContent(startPos:endPos),'format') /= 'binary') &
                call IO_error(error_ID = 844, ext_msg='format (coordinates)')
              dataType = getXMLValue(fileContent(startPos:endPos),'type')

              startPos = endPos + 2_pI64
              endPos  = startPos + index(fileContent(startPos:),IO_EOL,kind=pI64) - 2_pI64
              s = startPos + verify(fileContent(startPos:endPos),IO_WHITESPACE,kind=pI64) -1_pI64   ! start (no leading whitespace)

              coord = coord + 1

              call gridSizeOrigin(fileContent(s:endPos),headerType,compressed,dataType,coord)
            endif
            if(index(fileContent(startPos:endPos),'</Coordinates>',kind=pI64) /= 0_pI64) exit
            startPos = endPos + 2_pI64
          enddo
        endif
      endif
    endif

    if(gotCellData .and. gotCoordinates) exit
    startPos = endPos + 2_pI64

  end do
  material = material + 1
  if(.not. allocated(material))       call IO_error(error_ID = 844, ext_msg='material data not found')
  if(size(material) /= product(grid)) call IO_error(error_ID = 844, ext_msg='size(material)')
  if(any(geomSize<=0))                call IO_error(error_ID = 844, ext_msg='size')
  if(any(grid<1))                     call IO_error(error_ID = 844, ext_msg='grid')
  if(any(material<0))                 call IO_error(error_ID = 844, ext_msg='material ID < 0')

  contains

  !------------------------------------------------------------------------------------------------
  !> @brief determine size and origin from coordinates
  !------------------------------------------------------------------------------------------------
  subroutine gridSizeOrigin(base64_str,headerType,compressed,dataType,direction)

    character(len=*), intent(in) :: base64_str, &                                                   ! base64 encoded string of 1D coordinates
                                    headerType, &                                                   ! header type (UInt32 or Uint64)
                                    dataType                                                        ! data type (Int32, Int64, Float32, Float64)
    logical,          intent(in) :: compressed                                                      ! indicate whether data is zlib compressed
    integer,          intent(in) :: direction                                                       ! direction (1=x,2=y,3=z)

    real(pReal), dimension(:), allocatable :: coords,delta

    coords = as_pReal(base64_str,headerType,compressed,dataType)

    delta = coords(2:) - coords(:size(coords)-1)
    if(any(delta<0.0_pReal) .or. dNeq(maxval(delta),minval(delta),1.0e-8_pReal*maxval(abs(coords)))) &
      call IO_error(error_ID = 844, ext_msg = 'grid spacing')

    grid(direction)     = size(coords)-1
    origin(direction)   = coords(1)
    geomSize(direction) = coords(size(coords)) - coords(1)

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
        call IO_error(844_pInt,ext_msg='unknown data type: '//trim(dataType))
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
        call IO_error(844_pInt,ext_msg='unknown data type: '//trim(dataType))
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

    if(compressed) then
      bytes = asBytes_compressed(base64_str,headerType)
    else
      bytes = asBytes_uncompressed(base64_str,headerType)
    endif

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

    if    (headerType == 'UInt32') then
      temp = int(prec_bytesToC_INT32_T(base64_to_bytes(base64_str(:base64_nChar(4_pI64)))),pI64)
      nBlock = int(temp(1),pI64)
      headerLen = 4_pI64 * (3_pI64 + nBlock)
      temp = int(prec_bytesToC_INT32_T(base64_to_bytes(base64_str(:base64_nChar(headerLen)))),pI64)
    elseif(headerType == 'UInt64') then
      temp = int(prec_bytesToC_INT64_T(base64_to_bytes(base64_str(:base64_nChar(8_pI64)))),pI64)
      nBlock = int(temp(1),pI64)
      headerLen = 8_pI64 * (3_pI64 + nBlock)
      temp = int(prec_bytesToC_INT64_T(base64_to_bytes(base64_str(:base64_nChar(headerLen)))),pI64)
    endif

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
    enddo

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
    if    (headerType == 'UInt32') then
      do while(s+base64_nChar(4_pI64)<(len(base64_str,pI64)))
        nByte = int(prec_bytesToC_INT32_T(base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(4_pI64)))),pI64)
        bytes = [bytes,base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(4_pI64+nByte(1))),5_pI64)]
        s = s + base64_nChar(4_pI64+nByte(1))
      enddo
    elseif(headerType == 'UInt64') then
      do while(s+base64_nChar(8_pI64)<(len(base64_str,pI64)))
        nByte = int(prec_bytesToC_INT64_T(base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(8_pI64)))),pI64)
        bytes = [bytes,base64_to_bytes(base64_str(s+1_pI64:s+base64_nChar(8_pI64+nByte(1))),9_pI64)]
        s = s + base64_nChar(8_pI64+nByte(1))
      enddo
    endif

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
    if(s==0) then
      getXMLValue = ''
    else
      e = s + 1 + scan(line(s+1:),"'"//'"')
      if(scan(line(s:e-2),'=') == 0) then
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
      endif
    endif

  end function


  !------------------------------------------------------------------------------------------------
  !> @brief check for supported file format
  !------------------------------------------------------------------------------------------------
  pure function fileFormatOk(line)

    character(len=*),intent(in) :: line
    logical :: fileFormatOk

    fileFormatOk = getXMLValue(line,'type')       == 'RectilinearGrid' .and. &
                   getXMLValue(line,'byte_order') == 'LittleEndian' .and. &
                   getXMLValue(line,'compressor') /= 'vtkLZ4DataCompressor' .and. &
                   getXMLValue(line,'compressor') /= 'vtkLZMADataCompressor'

  end function fileFormatOk

end subroutine readVTR


!---------------------------------------------------------------------------------------------------
!> @brief Calculate undeformed position of IPs/cell centers (pretend to be an element)
!---------------------------------------------------------------------------------------------------
function IPcoordinates0(grid,geomSize,grid3Offset)

  integer,     dimension(3), intent(in) :: grid                                                     ! grid (for this process!)
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,                   intent(in) :: grid3Offset                                              ! grid(3) offset

  real(pReal), dimension(3,product(grid))  :: ipCoordinates0

  integer :: &
    a,b,c, &
    i

  i = 0
  do c = 1, grid(3); do b = 1, grid(2); do a = 1, grid(1)
    i = i + 1
    IPcoordinates0(1:3,i) = geomSize/real(grid,pReal) * (real([a,b,grid3Offset+c],pReal) -0.5_pReal)
  enddo; enddo; enddo

end function IPcoordinates0


!---------------------------------------------------------------------------------------------------
!> @brief Calculate position of undeformed nodes (pretend to be an element)
!---------------------------------------------------------------------------------------------------
pure function nodes0(grid,geomSize,grid3Offset)

  integer,     dimension(3), intent(in) :: grid                                                     ! grid (for this process!)
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,                   intent(in) :: grid3Offset                                              ! grid(3) offset

  real(pReal), dimension(3,product(grid+1)) :: nodes0

  integer :: &
    a,b,c, &
    n

  n = 0
  do c = 0, grid3; do b = 0, grid(2); do a = 0, grid(1)
    n = n + 1
    nodes0(1:3,n) = geomSize/real(grid,pReal) * real([a,b,grid3Offset+c],pReal)
  enddo; enddo; enddo

end function nodes0


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP interface areas
!--------------------------------------------------------------------------------------------------
pure function cellSurfaceArea(geomSize,grid)

  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,     dimension(3), intent(in) :: grid                                                     ! grid (for this process!)

  real(pReal), dimension(6,1,product(grid)) :: cellSurfaceArea

  cellSurfaceArea(1:2,1,:) = geomSize(2)/real(grid(2)) * geomSize(3)/real(grid(3))
  cellSurfaceArea(3:4,1,:) = geomSize(3)/real(grid(3)) * geomSize(1)/real(grid(1))
  cellSurfaceArea(5:6,1,:) = geomSize(1)/real(grid(1)) * geomSize(2)/real(grid(2))

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
pure function IPneighborhood(grid)

  integer, dimension(3), intent(in) :: grid                                                         ! grid (for this process!)

  integer, dimension(3,6,1,product(grid)) :: IPneighborhood                                         !< 6 neighboring IPs as [element ID, IP ID, face ID]

  integer :: &
   x,y,z, &
   e

  e = 0
  do z = 0,grid(3)-1; do y = 0,grid(2)-1; do x = 0,grid(1)-1
    e = e + 1
    ! element ID
    IPneighborhood(1,1,1,e) = z * grid(1) * grid(2) &
                            + y * grid(1) &
                            + modulo(x+1,grid(1)) &
                            + 1
    IPneighborhood(1,2,1,e) = z * grid(1) * grid(2) &
                            + y * grid(1) &
                            + modulo(x-1,grid(1)) &
                            + 1
    IPneighborhood(1,3,1,e) = z * grid(1) * grid(2) &
                            + modulo(y+1,grid(2)) * grid(1) &
                            + x &
                            + 1
    IPneighborhood(1,4,1,e) = z * grid(1) * grid(2) &
                            + modulo(y-1,grid(2)) * grid(1) &
                            + x &
                            + 1
    IPneighborhood(1,5,1,e) = modulo(z+1,grid(3)) * grid(1) * grid(2) &
                            + y * grid(1) &
                            + x &
                            + 1
    IPneighborhood(1,6,1,e) = modulo(z-1,grid(3)) * grid(1) * grid(2) &
                            + y * grid(1) &
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

  enddo; enddo; enddo

end function IPneighborhood


end module discretization_grid
