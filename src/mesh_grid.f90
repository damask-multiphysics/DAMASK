!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc, Abaqus and the spectral solver
!--------------------------------------------------------------------------------------------------
module mesh
  use, intrinsic :: iso_c_binding
  use prec
  use geometry_plastic_nonlocal
  use mesh_base
#include <petsc/finclude/petscsys.h>
  use PETScsys
  use DAMASK_interface
  use IO
  use debug
  use numerics
  use FEsolving

 
  implicit none
  private
 
  include 'fftw3-mpi.f03'
  integer, public, protected :: &
    mesh_Nnodes
 
  integer, dimension(:), allocatable, private :: &
    microGlobal
  integer, dimension(:), allocatable, private :: &
    mesh_homogenizationAt
 
  integer, dimension(:,:), allocatable, public, protected :: &
    mesh_element                                                                                    !< entryCount and list of elements containing node
 
  integer, dimension(:,:,:,:), allocatable, public, protected :: &
    mesh_ipNeighborhood                                                                             !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]
 
  real(pReal), public, protected :: &
    mesh_unitlength                                                                                 !< physical length of one unit in mesh
 
  real(pReal), dimension(:,:), allocatable, private :: &
    mesh_node                                                                                       !< node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)
 
 
  real(pReal), dimension(:,:), allocatable, public, protected :: &
    mesh_ipVolume, &                                                                                !< volume associated with IP (initially!)
    mesh_node0                                                                                      !< node x,y,z coordinates (initially!)
 
  real(pReal), dimension(:,:,:), allocatable, public, protected :: &
    mesh_ipArea                                                                                     !< area of interface to neighboring IP (initially!)
 
  real(pReal), dimension(:,:,:), allocatable, public :: &
    mesh_ipCoordinates                                                                              !< IP x,y,z coordinates (after deformation!)
 
  real(pReal),dimension(:,:,:,:), allocatable, public, protected :: &
    mesh_ipAreaNormal                                                                               !< area normal of interface to neighboring IP (initially!)
 
  logical, dimension(3), public, parameter :: mesh_periodicSurface = .true.                         !< flag indicating periodic outer surfaces (used for fluxes)
 
 
 ! grid specific
  integer, dimension(3), public, protected :: &
    grid                                                                                            !< (global) grid
  integer, public, protected :: &
    mesh_NcpElemsGlobal, &                                                                          !< total number of CP elements in global mesh
    grid3, &                                                                                        !< (local) grid in 3rd direction
    grid3Offset                                                                                     !< (local) grid offset in 3rd direction
  real(pReal), dimension(3), public, protected :: &
    geomSize
  real(pReal), public, protected :: &
    size3, &                                                                                        !< (local) size in 3rd direction
    size3offset                                                                                     !< (local) size offset in 3rd direction
 
  public :: &
    mesh_init
 
  private :: &
    mesh_build_ipAreas, &
    mesh_build_ipNormals, &
    mesh_spectral_build_nodes, &
    mesh_spectral_build_elements, &
    mesh_spectral_build_ipNeighborhood, &
    mesh_build_ipCoordinates
 
  type, public, extends(tMesh) :: tMesh_grid
  
   integer, dimension(3), public :: &
    grid                                                                                             !< (global) grid
  integer, public :: &
    mesh_NcpElemsGlobal, &                                                                           !< total number of CP elements in global mesh
    grid3, &                                                                                         !< (local) grid in 3rd direction
    grid3Offset                                                                                      !< (local) grid offset in 3rd direction
  real(pReal), dimension(3), public :: &
    geomSize
  real(pReal), public :: &
    size3, &                                                                                         !< (local) size in 3rd direction
    size3offset
    
    contains
    procedure, pass(self) :: tMesh_grid_init
    generic, public :: init => tMesh_grid_init
  end type tMesh_grid
  
  type(tMesh_grid), public, protected :: theMesh
 
contains

subroutine tMesh_grid_init(self,nodes)
 
 class(tMesh_grid) :: self
 real(pReal), dimension(:,:), intent(in) :: nodes
 
 call self%tMesh%init('grid',10,nodes)
 
end subroutine tMesh_grid_init

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

  integer(C_INTPTR_T) :: devNull, local_K, local_K_offset
  integer :: ierr, worldsize, j
  integer, intent(in), optional :: el, ip
  logical :: myDebug

  write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'

  mesh_unitlength = numerics_unitlength                                                             ! set physical extent of a length unit in mesh

  myDebug = (iand(debug_level(debug_mesh),debug_levelBasic) /= 0)

  call fftw_mpi_init()
  call mesh_spectral_read_grid()


  call MPI_comm_size(PETSC_COMM_WORLD, worldsize, ierr)
  if(ierr /=0) call IO_error(894, ext_msg='MPI_comm_size')
  if(worldsize>grid(3)) call IO_error(894, ext_msg='number of processes exceeds grid(3)')


  devNull = fftw_mpi_local_size_3d(int(grid(3),C_INTPTR_T), &
                                   int(grid(2),C_INTPTR_T), &
                                   int(grid(1),C_INTPTR_T)/2+1, &
                                   PETSC_COMM_WORLD, &
                                   local_K, &                                                       ! domain grid size along z
                                   local_K_offset)                                                  ! domain grid offset along z
  grid3       = int(local_K,pInt)
  grid3Offset = int(local_K_offset,pInt)
  size3       = geomSize(3)*real(grid3,pReal)      /real(grid(3),pReal)
  size3Offset = geomSize(3)*real(grid3Offset,pReal)/real(grid(3),pReal)

  mesh_NcpElemsGlobal = product(grid)

  mesh_Nnodes  = product(grid(1:2) + 1)*(grid3 + 1)

  mesh_node0 = mesh_spectral_build_nodes()
  mesh_node  = mesh_node0
  if (myDebug) write(6,'(a)') ' Built nodes'; flush(6)

  call theMesh%init(mesh_node)
  call theMesh%setNelems(product(grid(1:2))*grid3)
  call mesh_spectral_build_elements()
  mesh_homogenizationAt = mesh_homogenizationAt(product(grid(1:2))*grid3Offset+1: &
                                                product(grid(1:2))*(grid3Offset+grid3))             ! reallocate/shrink in case of MPI
  
  if (myDebug) write(6,'(a)') ' Built elements'; flush(6)
  
  

  if (myDebug) write(6,'(a)') ' Built cell nodes'; flush(6)
  mesh_ipCoordinates = mesh_build_ipCoordinates()
  if (myDebug) write(6,'(a)') ' Built IP coordinates'; flush(6)
  allocate(mesh_ipVolume(1,theMesh%nElems),source=product([geomSize(1:2),size3]/real([grid(1:2),grid3])))
  if (myDebug) write(6,'(a)') ' Built IP volumes'; flush(6)
  mesh_ipArea       = mesh_build_ipAreas()
  mesh_ipAreaNormal = mesh_build_ipNormals()
  if (myDebug) write(6,'(a)') ' Built IP areas'; flush(6)

  call mesh_spectral_build_ipNeighborhood
  call geometry_plastic_nonlocal_set_IPneighborhood(mesh_ipNeighborhood)

  if (myDebug) write(6,'(a)') ' Built IP neighborhood'; flush(6)

  if (debug_e < 1 .or. debug_e > theMesh%nElems) &
    call IO_error(602,ext_msg='element')                                                            ! selected element does not exist
  if (debug_i < 1 .or. debug_i > theMesh%elem%nIPs) &
    call IO_error(602,ext_msg='IP')                                                                 ! selected element does not have requested IP

  FEsolving_execElem = [ 1,theMesh%nElems ]                                                         ! parallel loop bounds set to comprise all DAMASK elements
  allocate(FEsolving_execIP(2,theMesh%nElems), source=1)                                            ! parallel loop bounds set to comprise from first IP...
  forall (j = 1:theMesh%nElems) FEsolving_execIP(2,j) = theMesh%elem%nIPs                           ! ...up to own IP count for each element


!!!! COMPATIBILITY HACK !!!!
  theMesh%homogenizationAt  = mesh_element(3,:)
  theMesh%microstructureAt  = mesh_element(4,:)
!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Parses geometry file
!> @details important variables have an implicit "save" attribute. Therefore, this function is 
! supposed to be called only once!
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_read_grid()

  character(len=:),            allocatable :: rawData
  character(len=65536)                     :: line
  integer, allocatable, dimension(:) :: chunkPos
  integer :: h =- 1
  integer ::  &
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
! read data as stream
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
    call IO_error(error_ID=841, ext_msg='mesh_spectral_read_grid')
  else
    chunkPos = IO_stringPos(rawData(1:endPos))
    if (chunkPos(1) < 2) call IO_error(error_ID=841, ext_msg='mesh_spectral_read_grid')
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
    call IO_error(error_ID = 842, ext_msg='homogenization (mesh_spectral_read_grid)')
  if(any(grid < 1)) &
    call IO_error(error_ID = 842, ext_msg='grid (mesh_spectral_read_grid)')
  if(any(geomSize < 0.0_pReal)) &
    call IO_error(error_ID = 842, ext_msg='size (mesh_spectral_read_grid)')

  allocate(microGlobal(product(grid)), source = -1)
  allocate(mesh_homogenizationAt(product(grid)), source = h)                                        ! too large in case of MPI (shrink later, not very elegant)
     
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
      microGlobal(e:e+c-1) =  [(IO_intValue(line,chunkPos,i+1), i=0, c-1)]
    else noCompression
      compression: if (IO_lc(IO_stringValue(line,chunkPos,2))  == 'of') then
        c = IO_intValue(line,chunkPos,1)
        microGlobal(e:e+c-1) = [(IO_intValue(line,chunkPos,3),i = 1,IO_intValue(line,chunkPos,1))]
      else if (IO_lc(IO_stringValue(line,chunkPos,2))  == 'to') then compression
        c = abs(IO_intValue(line,chunkPos,3) - IO_intValue(line,chunkPos,1)) + 1
        o = merge(+1, -1, IO_intValue(line,chunkPos,3) > IO_intValue(line,chunkPos,1))
        microGlobal(e:e+c-1) = [(i, i = IO_intValue(line,chunkPos,1),IO_intValue(line,chunkPos,3),o)]
      else compression
        c = chunkPos(1)
        microGlobal(e:e+c-1) =  [(IO_intValue(line,chunkPos,i+1), i=0, c-1)]
      endif compression
    endif noCompression

    e = e+c
  end do

  if (e-1 /= product(grid)) call IO_error(error_ID = 843, el=e)

end subroutine mesh_spectral_read_grid


!---------------------------------------------------------------------------------------------------
!> @brief Calculates position of nodes (pretend to be an element)
!---------------------------------------------------------------------------------------------------
pure function mesh_spectral_build_nodes()

  real(pReal), dimension(3,mesh_Nnodes) :: mesh_spectral_build_nodes
  integer :: n,a,b,c

  n = 0
  do c = 0, grid3
    do b = 0, grid(2)
      do a = 0, grid(1)
         n = n + 1
         mesh_spectral_build_nodes(1:3,n) = geomSize/real(grid,pReal) * real([a,b,grid3Offset+c],pReal)
      enddo
    enddo
  enddo

end function mesh_spectral_build_nodes


!---------------------------------------------------------------------------------------------------
!> @brief Calculates position of IPs/cell centres (pretend to be an element)
!---------------------------------------------------------------------------------------------------
function mesh_build_ipCoordinates()

  real(pReal), dimension(3,1,theMesh%nElems) :: mesh_build_ipCoordinates
  integer :: n,a,b,c
 
  n = 0
  do c = 1, grid3
    do b = 1, grid(2)
      do a = 1, grid(1)
         n = n + 1
         mesh_build_ipCoordinates(1:3,1,n) = geomSize/real(grid,pReal) * (real([a,b,grid3Offset+c],pReal) -0.5_pReal)
      enddo
    enddo
  enddo

end function mesh_build_ipCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, material, texture, and node list per element.
!! Allocates global array 'mesh_element'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_elements

  integer :: &
    e, &
    elemOffset

  allocate(mesh_element    (4+8,theMesh%nElems), source = 0)

  elemOffset = product(grid(1:2))*grid3Offset
  do e=1, theMesh%nElems
    mesh_element( 1,e) = -1                                                                         ! DEPRECATED
    mesh_element( 2,e) = -1                                                                         ! DEPRECATED
    mesh_element( 3,e) = mesh_homogenizationAt(e)
    mesh_element( 4,e) = microGlobal(e+elemOffset)                                                  ! microstructure
    mesh_element( 5,e) = e + (e-1)/grid(1) + &
                                      ((e-1)/(grid(1)*grid(2)))*(grid(1)+1)                         ! base node
    mesh_element( 6,e) = mesh_element(5,e) + 1
    mesh_element( 7,e) = mesh_element(5,e) + grid(1) + 2
    mesh_element( 8,e) = mesh_element(5,e) + grid(1) + 1
    mesh_element( 9,e) = mesh_element(5,e) +(grid(1) + 1) * (grid(2) + 1)                           ! second floor base node
    mesh_element(10,e) = mesh_element(9,e) + 1
    mesh_element(11,e) = mesh_element(9,e) + grid(1) + 2
    mesh_element(12,e) = mesh_element(9,e) + grid(1) + 1
  enddo

end subroutine mesh_spectral_build_elements


!--------------------------------------------------------------------------------------------------
!> @brief build neighborhood relations for spectral
!> @details assign globals: mesh_ipNeighborhood
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_ipNeighborhood

  integer :: &
   x,y,z, &
   e
  allocate(mesh_ipNeighborhood(3,6,1,theMesh%nElems),source=0)

  e = 0
  do z = 0,grid3-1
    do y = 0,grid(2)-1
      do x = 0,grid(1)-1
        e = e + 1
          ! neigboring element
          mesh_ipNeighborhood(1,1,1,e) = z * grid(1) * grid(2) &
                                       + y * grid(1) &
                                       + modulo(x+1,grid(1)) &
                                       + 1
          mesh_ipNeighborhood(1,2,1,e) = z * grid(1) * grid(2) &
                                       + y * grid(1) &
                                       + modulo(x-1,grid(1)) &
                                       + 1
          mesh_ipNeighborhood(1,3,1,e) = z * grid(1) * grid(2) &
                                       + modulo(y+1,grid(2)) * grid(1) &
                                       + x &
                                       + 1
          mesh_ipNeighborhood(1,4,1,e) = z * grid(1) * grid(2) &
                                       + modulo(y-1,grid(2)) * grid(1) &
                                       + x &
                                       + 1
          mesh_ipNeighborhood(1,5,1,e) = modulo(z+1,grid3) * grid(1) * grid(2) &
                                       + y * grid(1) &
                                       + x &
                                       + 1
          mesh_ipNeighborhood(1,6,1,e) = modulo(z-1,grid3) * grid(1) * grid(2) &
                                       + y * grid(1) &
                                       + x &
                                       + 1
          ! neigboring IP
          mesh_ipNeighborhood(2,1:6,1,e) = 1
          ! neigboring face
          mesh_ipNeighborhood(3,1,1,e) = 2
          mesh_ipNeighborhood(3,2,1,e) = 1
          mesh_ipNeighborhood(3,3,1,e) = 4
          mesh_ipNeighborhood(3,4,1,e) = 3
          mesh_ipNeighborhood(3,5,1,e) = 6
          mesh_ipNeighborhood(3,6,1,e) = 5
      enddo
    enddo
  enddo

end subroutine mesh_spectral_build_ipNeighborhood


!--------------------------------------------------------------------------------------------------
!> @brief builds mesh of (distorted) cubes for given coordinates (= center of the cubes)
!--------------------------------------------------------------------------------------------------
function mesh_nodesAroundCentres(gDim,Favg,centres) result(nodes)
 
  real(pReal), intent(in), dimension(:,:,:,:) :: &
    centres
  real(pReal),             dimension(3,size(centres,2)+1,size(centres,3)+1,size(centres,4)+1) :: &
    nodes
  real(pReal), intent(in), dimension(3) :: &
    gDim
  real(pReal), intent(in), dimension(3,3) :: &
    Favg
  real(pReal),             dimension(3,size(centres,2)+2,size(centres,3)+2,size(centres,4)+2) :: &
    wrappedCentres
 
  integer :: &
    i,j,k,n
  integer,           dimension(3), parameter :: &
    diag = 1
  integer,           dimension(3) :: &
    shift = 0, &
    lookup = 0, &
    me = 0, &
    iRes = 0
  integer,           dimension(3,8) :: &
    neighbor = reshape([ &
                        0, 0, 0, &
                        1, 0, 0, &
                        1, 1, 0, &
                        0, 1, 0, &
                        0, 0, 1, &
                        1, 0, 1, &
                        1, 1, 1, &
                        0, 1, 1  ], [3,8])

!--------------------------------------------------------------------------------------------------
! initializing variables
 iRes =  [size(centres,2),size(centres,3),size(centres,4)]
 nodes = 0.0_pReal
 wrappedCentres = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! building wrappedCentres = centroids + ghosts
  wrappedCentres(1:3,2:iRes(1)+1,2:iRes(2)+1,2:iRes(3)+1) = centres
  do k = 0,iRes(3)+1
    do j = 0,iRes(2)+1
      do i = 0,iRes(1)+1
        if (k==0 .or. k==iRes(3)+1 .or. &                                                           ! z skin
            j==0 .or. j==iRes(2)+1 .or. &                                                           ! y skin
            i==0 .or. i==iRes(1)+1      ) then                                                      ! x skin
          me = [i,j,k]                                                                              ! me on skin
          shift = sign(abs(iRes+diag-2*me)/(iRes+diag),iRes+diag-2*me)
          lookup = me-diag+shift*iRes
          wrappedCentres(1:3,i+1,        j+1,        k+1) = &
                 centres(1:3,lookup(1)+1,lookup(2)+1,lookup(3)+1) &
                 - matmul(Favg, real(shift,pReal)*gDim)
        endif
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! averaging
  do k = 0,iRes(3); do j = 0,iRes(2); do i = 0,iRes(1)
    do n = 1,8
     nodes(1:3,i+1,j+1,k+1) = &
     nodes(1:3,i+1,j+1,k+1) + wrappedCentres(1:3,i+1+neighbor(1,n), &
                                                                j+1+neighbor(2,n), &
                                                                k+1+neighbor(3,n) )
    enddo
  enddo; enddo; enddo
  nodes = nodes/8.0_pReal

end function mesh_nodesAroundCentres


!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
pure function mesh_build_ipAreas()

  real(pReal), dimension(6,1,theMesh%nElems) :: mesh_build_ipAreas

  mesh_build_ipAreas(1:2,1,:) = geomSize(2)/real(grid(2)) * geomSize(3)/real(grid(3))
  mesh_build_ipAreas(3:4,1,:) = geomSize(3)/real(grid(3)) * geomSize(1)/real(grid(1))
  mesh_build_ipAreas(5:6,1,:) = geomSize(1)/real(grid(1)) * geomSize(2)/real(grid(2))
  
end function mesh_build_ipAreas


!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
pure function mesh_build_ipNormals()

  real, dimension(3,6,1,theMesh%nElems) :: mesh_build_ipNormals

  mesh_build_ipNormals(1:3,1,1,:) = spread([+1.0_pReal, 0.0_pReal, 0.0_pReal],2,theMesh%nElems)
  mesh_build_ipNormals(1:3,2,1,:) = spread([-1.0_pReal, 0.0_pReal, 0.0_pReal],2,theMesh%nElems)
  mesh_build_ipNormals(1:3,3,1,:) = spread([ 0.0_pReal,+1.0_pReal, 0.0_pReal],2,theMesh%nElems)
  mesh_build_ipNormals(1:3,4,1,:) = spread([ 0.0_pReal,-1.0_pReal, 0.0_pReal],2,theMesh%nElems)
  mesh_build_ipNormals(1:3,5,1,:) = spread([ 0.0_pReal, 0.0_pReal,+1.0_pReal],2,theMesh%nElems)
  mesh_build_ipNormals(1:3,6,1,:) = spread([ 0.0_pReal, 0.0_pReal,-1.0_pReal],2,theMesh%nElems)
  
end function mesh_build_ipNormals


end module mesh
