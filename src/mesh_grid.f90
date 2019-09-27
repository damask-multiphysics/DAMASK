!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse geometry file to set up discretization and geometry for nonlocal model
!--------------------------------------------------------------------------------------------------
module mesh_grid
#include <petsc/finclude/petscsys.h>
  use PETScsys

  use prec
  use system_routines
  use DAMASK_interface
  use IO
  use debug
  use numerics
  use discretization
  use geometry_plastic_nonlocal
  use FEsolving
 
  implicit none
  private
  
  integer,     dimension(3), public, protected :: &
    grid                                                                                            !< (global) grid
  integer,                   public, protected :: &
    grid3, &                                                                                        !< (local) grid in 3rd direction
    grid3Offset                                                                                     !< (local) grid offset in 3rd direction
    
  real(pReal), dimension(3), public, protected :: &
    geomSize
  real(pReal),               public, protected :: &
    size3, &                                                                                        !< (local) size in 3rd direction
    size3offset                                                                                     !< (local) size offset in 3rd direction
 
  public :: &
    mesh_init
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief reads the geometry file to obtain information on discretization
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

  integer, intent(in), optional :: el, ip                                                           ! for compatibility reasons
  
  include 'fftw3-mpi.f03'
  real(pReal), dimension(3) :: &
    mySize                                                                                          !< domain size of this process
  integer,     dimension(3) :: &
    myGrid                                                                                          !< domain grid of this process

  integer,     dimension(:),   allocatable :: &
    microstructureAt, &
    homogenizationAt

  integer :: j
  integer(C_INTPTR_T) :: &
    devNull, z, z_offset

  write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'

  call readGeom(grid,geomSize,microstructureAt,homogenizationAt)

!--------------------------------------------------------------------------------------------------
! grid solver specific quantities
  if(worldsize>grid(3)) call IO_error(894, ext_msg='number of processes exceeds grid(3)')

  call fftw_mpi_init
  devNull = fftw_mpi_local_size_3d(int(grid(3),C_INTPTR_T), &
                                   int(grid(2),C_INTPTR_T), &
                                   int(grid(1),C_INTPTR_T)/2+1, &
                                   PETSC_COMM_WORLD, &
                                   z, &                                                             ! domain grid size along z
                                   z_offset)                                                        ! domain grid offset along z
  grid3       = int(z)
  grid3Offset = int(z_offset)
  size3       = geomSize(3)*real(grid3,pReal)      /real(grid(3),pReal)
  size3Offset = geomSize(3)*real(grid3Offset,pReal)/real(grid(3),pReal)
  myGrid = [grid(1:2),grid3]
  mySize = [geomSize(1:2),size3]

!--------------------------------------------------------------------------------------------------
! general discretization
  microstructureAt = microstructureAt(product(grid(1:2))*grid3Offset+1: &
                                      product(grid(1:2))*(grid3Offset+grid3))                       ! reallocate/shrink in case of MPI
  homogenizationAt = homogenizationAt(product(grid(1:2))*grid3Offset+1: &
                                      product(grid(1:2))*(grid3Offset+grid3))                       ! reallocate/shrink in case of MPI

  call discretization_init(homogenizationAt,microstructureAt, &
                           reshape(IPcoordinates(myGrid,mySize,grid3Offset),[3,product(myGrid)]), &
                           Nodes(myGrid,mySize,grid3Offset))

  FEsolving_execElem = [1,product(myGrid)]                                                          ! parallel loop bounds set to comprise all elements
  allocate(FEsolving_execIP(2,product(myGrid)),source=1)                                            ! parallel loop bounds set to comprise the only IP

!--------------------------------------------------------------------------------------------------
! geometry information required by the nonlocal CP model
  call geometry_plastic_nonlocal_setIPvolume(reshape([(product(mySize/real(myGrid,pReal)),j=1,product(myGrid))], &
                                                     [1,product(myGrid)]))
  call geometry_plastic_nonlocal_setIParea        (cellEdgeArea(mySize,myGrid))
  call geometry_plastic_nonlocal_setIPareaNormal  (cellEdgeNormal(product(myGrid)))
  call geometry_plastic_nonlocal_setIPneighborhood(IPneighborhood(myGrid))

!--------------------------------------------------------------------------------------------------
! sanity checks for debugging
  if (debug_e < 1 .or. debug_e > product(myGrid)) call IO_error(602,ext_msg='element')              ! selected element does not exist
  if (debug_i /= 1)                               call IO_error(602,ext_msg='IP')                   ! selected IP does not exist

end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Parses geometry file
!> @details important variables have an implicit "save" attribute. Therefore, this function is 
! supposed to be called only once!
!--------------------------------------------------------------------------------------------------
subroutine readGeom(grid,geomSize,microstructure,homogenization)

  integer,     dimension(3), intent(out)              :: grid                                       ! grid (for all processes!)
  real(pReal), dimension(3), intent(out)              :: geomSize                                   ! size (for all processes!)
  integer,     dimension(:), intent(out), allocatable :: &
    microstructure, &
    homogenization
   
  character(len=:),      allocatable :: rawData
  character(len=65536)               :: line
  integer, allocatable, dimension(:) :: chunkPos
  integer :: &
    h =- 1, &
    headerLength = -1, &                                                                            !< length of header (in lines)
    fileLength, &                                                                                   !< length of the geom file (in characters)
    fileUnit, &
    startPos, endPos, &
    myStat, &
    l, &                                                                                            !< line counter
    c, &                                                                                            !< counter for # microstructures in line
    o, &                                                                                            !< order of "to" packing
    e, &                                                                                            !< "element", i.e. spectral collocation point 
    i, j
    
  grid = -1
  geomSize = -1.0_pReal

!--------------------------------------------------------------------------------------------------
! read raw data as stream
  inquire(file = trim(geometryFile), size=fileLength)
  open(newunit=fileUnit, file=trim(geometryFile), access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if(myStat /= 0) call IO_error(100,ext_msg=trim(geometryFile))
  allocate(character(len=fileLength)::rawData)
  read(fileUnit) rawData
  close(fileUnit)
  
!--------------------------------------------------------------------------------------------------
! get header length
  endPos = index(rawData,new_line(''))
  if(endPos <= index(rawData,'head')) then
    startPos = len(rawData)
    call IO_error(error_ID=841, ext_msg='readGeom')
  else
    chunkPos = IO_stringPos(rawData(1:endPos))
    if (chunkPos(1) < 2) call IO_error(error_ID=841, ext_msg='readGeom')
    headerLength = IO_intValue(rawData(1:endPos),chunkPos,1)
    startPos = endPos + 1
  endif

!--------------------------------------------------------------------------------------------------
! read and interprete header
  l = 0
  do while (l < headerLength .and. startPos < len(rawData))
    endPos = startPos + index(rawData(startPos:),new_line('')) - 1
    if (endPos < startPos) endPos = len(rawData)                                                    ! end of file without new line
    line = rawData(startPos:endPos)
    startPos = endPos + 1
    l = l + 1

    chunkPos = IO_stringPos(trim(line))
    if (chunkPos(1) < 2) cycle                                                                      ! need at least one keyword value pair
    
    select case ( IO_lc(IO_StringValue(trim(line),chunkPos,1,.true.)) )
      case ('grid')
        if (chunkPos(1) > 6) then
          do j = 2,6,2
            select case (IO_lc(IO_stringValue(line,chunkPos,j)))
              case('a')
                grid(1) = IO_intValue(line,chunkPos,j+1)
              case('b')
                grid(2) = IO_intValue(line,chunkPos,j+1)
              case('c')
                grid(3) = IO_intValue(line,chunkPos,j+1)
            end select
          enddo
        endif
        
      case ('size')
        if (chunkPos(1) > 6) then
          do j = 2,6,2
            select case (IO_lc(IO_stringValue(line,chunkPos,j)))
              case('x')
                geomSize(1) = IO_floatValue(line,chunkPos,j+1)
              case('y')
                geomSize(2) = IO_floatValue(line,chunkPos,j+1)
              case('z')
                geomSize(3) = IO_floatValue(line,chunkPos,j+1)
            end select
          enddo
        endif
        
      case ('homogenization')
        if (chunkPos(1) > 1) h = IO_intValue(line,chunkPos,2)
    end select

  enddo

!--------------------------------------------------------------------------------------------------
! sanity checks
  if(h < 1) &
    call IO_error(error_ID = 842, ext_msg='homogenization (readGeom)')
  if(any(grid < 1)) &
    call IO_error(error_ID = 842, ext_msg='grid (readGeom)')
  if(any(geomSize < 0.0_pReal)) &
    call IO_error(error_ID = 842, ext_msg='size (readGeom)')

  allocate(microstructure(product(grid)), source = -1)                                              ! too large in case of MPI (shrink later, not very elegant)
  allocate(homogenization(product(grid)), source = h)                                               ! too large in case of MPI (shrink later, not very elegant)
     
!--------------------------------------------------------------------------------------------------
! read and interpret content
  e = 1
  do while (startPos < len(rawData))
    endPos = startPos + index(rawData(startPos:),new_line('')) - 1
    if (endPos < startPos) endPos = len(rawData)                                                    ! end of file without new line
    line = rawData(startPos:endPos)
    startPos = endPos + 1
    l = l + 1
    chunkPos = IO_stringPos(trim(line))
    
    noCompression: if (chunkPos(1) /= 3) then
      c = chunkPos(1)
      microstructure(e:e+c-1) =  [(IO_intValue(line,chunkPos,i+1), i=0, c-1)]
    else noCompression
      compression: if (IO_lc(IO_stringValue(line,chunkPos,2))  == 'of') then
        c = IO_intValue(line,chunkPos,1)
        microstructure(e:e+c-1) = [(IO_intValue(line,chunkPos,3),i = 1,IO_intValue(line,chunkPos,1))]
      else if (IO_lc(IO_stringValue(line,chunkPos,2))  == 'to') then compression
        c = abs(IO_intValue(line,chunkPos,3) - IO_intValue(line,chunkPos,1)) + 1
        o = merge(+1, -1, IO_intValue(line,chunkPos,3) > IO_intValue(line,chunkPos,1))
        microstructure(e:e+c-1) = [(i, i = IO_intValue(line,chunkPos,1),IO_intValue(line,chunkPos,3),o)]
      else compression
        c = chunkPos(1)
        microstructure(e:e+c-1) =  [(IO_intValue(line,chunkPos,i+1), i=0, c-1)]
      endif compression
    endif noCompression

    e = e+c
  end do

  if (e-1 /= product(grid)) call IO_error(error_ID = 843, el=e)

end subroutine readGeom


!---------------------------------------------------------------------------------------------------
!> @brief Calculate position of IPs/cell centres (pretend to be an element)
!---------------------------------------------------------------------------------------------------
function IPcoordinates(grid,geomSize,grid3Offset)

  integer,     dimension(3), intent(in) :: grid                                                     ! grid (for this process!)
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,                   intent(in) :: grid3Offset                                              ! grid(3) offset

  real(pReal), dimension(3,1,product(grid))  :: ipCoordinates

  integer :: &
    a,b,c, &
    i
 
  i = 0
  do c = 1, grid(3); do b = 1, grid(2); do a = 1, grid(1)
    i = i + 1
    IPcoordinates(1:3,1,i) = geomSize/real(grid,pReal) * (real([a,b,grid3Offset+c],pReal) -0.5_pReal)
  enddo; enddo; enddo

end function IPcoordinates


!---------------------------------------------------------------------------------------------------
!> @brief Calculate position of nodes (pretend to be an element)
!---------------------------------------------------------------------------------------------------
pure function nodes(grid,geomSize,grid3Offset)

  integer,     dimension(3), intent(in) :: grid                                                     ! grid (for this process!)
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,                   intent(in) :: grid3Offset                                              ! grid(3) offset

  real(pReal), dimension(3,product(grid+1))  :: nodes

  integer :: &
    a,b,c, &
    n

  n = 0
  do c = 0, grid3; do b = 0, grid(2); do a = 0, grid(1)
    n = n + 1
    nodes(1:3,n) = geomSize/real(grid,pReal) * real([a,b,grid3Offset+c],pReal)
  enddo; enddo; enddo

end function nodes


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP interface areas
!--------------------------------------------------------------------------------------------------
pure function cellEdgeArea(geomSize,grid)
  
  real(pReal), dimension(3), intent(in) :: geomSize                                                 ! size (for this process!)
  integer,     dimension(3), intent(in) :: grid                                                     ! grid (for this process!)
  
  real(pReal), dimension(6,1,product(grid)) :: cellEdgeArea

  cellEdgeArea(1:2,1,:) = geomSize(2)/real(grid(2)) * geomSize(3)/real(grid(3))
  cellEdgeArea(3:4,1,:) = geomSize(3)/real(grid(3)) * geomSize(1)/real(grid(1))
  cellEdgeArea(5:6,1,:) = geomSize(1)/real(grid(1)) * geomSize(2)/real(grid(2))
  
end function cellEdgeArea


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP interface areas normals
!--------------------------------------------------------------------------------------------------
pure function cellEdgeNormal(nElems)

  integer, intent(in) :: nElems
  
  real, dimension(3,6,1,nElems) :: cellEdgeNormal

  cellEdgeNormal(1:3,1,1,:) = spread([+1.0_pReal, 0.0_pReal, 0.0_pReal],2,nElems)
  cellEdgeNormal(1:3,2,1,:) = spread([-1.0_pReal, 0.0_pReal, 0.0_pReal],2,nElems)
  cellEdgeNormal(1:3,3,1,:) = spread([ 0.0_pReal,+1.0_pReal, 0.0_pReal],2,nElems)
  cellEdgeNormal(1:3,4,1,:) = spread([ 0.0_pReal,-1.0_pReal, 0.0_pReal],2,nElems)
  cellEdgeNormal(1:3,5,1,:) = spread([ 0.0_pReal, 0.0_pReal,+1.0_pReal],2,nElems)
  cellEdgeNormal(1:3,6,1,:) = spread([ 0.0_pReal, 0.0_pReal,-1.0_pReal],2,nElems)
  
end function cellEdgeNormal


!--------------------------------------------------------------------------------------------------
!> @brief Build IP neighborhood relations
!--------------------------------------------------------------------------------------------------
pure function IPneighborhood(grid)

  integer, dimension(3), intent(in) :: grid                                                         ! grid (for this process!)
  
  integer, dimension(3,6,1,product(grid)) :: IPneighborhood                                         !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]
  
  integer :: &
   x,y,z, &
   e

  e = 0
  do z = 0,grid(3)-1; do y = 0,grid(2)-1; do x = 0,grid(1)-1
    e = e + 1
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
    IPneighborhood(2,1:6,1,e) = 1
    IPneighborhood(3,1,1,e) = 2
    IPneighborhood(3,2,1,e) = 1
    IPneighborhood(3,3,1,e) = 4
    IPneighborhood(3,4,1,e) = 3
    IPneighborhood(3,5,1,e) = 6
    IPneighborhood(3,6,1,e) = 5
  enddo; enddo; enddo
 
end function IPneighborhood


!!--------------------------------------------------------------------------------------------------
!!> @brief builds mesh of (distorted) cubes for given coordinates (= center of the cubes)
!!--------------------------------------------------------------------------------------------------
!function mesh_nodesAroundCentres(gDim,Favg,centres) result(nodes)
!
!  real(pReal), intent(in), dimension(:,:,:,:) :: &
!    centres
!  real(pReal),             dimension(3,size(centres,2)+1,size(centres,3)+1,size(centres,4)+1) :: &
!    nodes
!  real(pReal), intent(in), dimension(3) :: &
!    gDim
!  real(pReal), intent(in), dimension(3,3) :: &
!    Favg
!  real(pReal),             dimension(3,size(centres,2)+2,size(centres,3)+2,size(centres,4)+2) :: &
!    wrappedCentres
!
!  integer :: &
!    i,j,k,n
!  integer,           dimension(3), parameter :: &
!    diag = 1
!  integer,           dimension(3) :: &
!    shift = 0, &
!    lookup = 0, &
!    me = 0, &
!    iRes = 0
!  integer,           dimension(3,8) :: &
!    neighbor = reshape([ &
!                        0, 0, 0, &
!                        1, 0, 0, &
!                        1, 1, 0, &
!                        0, 1, 0, &
!                        0, 0, 1, &
!                        1, 0, 1, &
!                        1, 1, 1, &
!                        0, 1, 1  ], [3,8])
!
!!--------------------------------------------------------------------------------------------------
!! initializing variables
! iRes =  [size(centres,2),size(centres,3),size(centres,4)]
! nodes = 0.0_pReal
! wrappedCentres = 0.0_pReal
!
!!--------------------------------------------------------------------------------------------------
!! building wrappedCentres = centroids + ghosts
!  wrappedCentres(1:3,2:iRes(1)+1,2:iRes(2)+1,2:iRes(3)+1) = centres
!  do k = 0,iRes(3)+1
!    do j = 0,iRes(2)+1
!      do i = 0,iRes(1)+1
!        if (k==0 .or. k==iRes(3)+1 .or. &                                                  ! z skin
!            j==0 .or. j==iRes(2)+1 .or. &                                                  ! y skin
!            i==0 .or. i==iRes(1)+1      ) then                                             ! x skin
!          me = [i,j,k]                                                                               ! me on skin
!          shift = sign(abs(iRes+diag-2*me)/(iRes+diag),iRes+diag-2*me)
!          lookup = me-diag+shift*iRes
!          wrappedCentres(1:3,i+1,        j+1,        k+1) = &
!                 centres(1:3,lookup(1)+1,lookup(2)+1,lookup(3)+1) &
!                 - matmul(Favg, real(shift,pReal)*gDim)
!        endif
!  enddo; enddo; enddo
!
!!--------------------------------------------------------------------------------------------------
!! averaging
!  do k = 0,iRes(3); do j = 0,iRes(2); do i = 0,iRes(1)
!    do n = 1,8
!     nodes(1:3,i+1,j+1,k+1) = &
!     nodes(1:3,i+1,j+1,k+1) + wrappedCentres(1:3,i+1+neighbor(1,n), &
!                                                                j+1+neighbor(2,n), &
!                                                                k+1+neighbor(3,n) )
!    enddo
!  enddo; enddo; enddo
!  nodes = nodes/8.0_pReal
!
!end function mesh_nodesAroundCentres

end module mesh_grid
