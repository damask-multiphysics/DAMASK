!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc, Abaqus and the spectral solver
!--------------------------------------------------------------------------------------------------
module mesh
#include <petsc/finclude/petscsys.h>
 use, intrinsic :: iso_c_binding
 use prec
 use debug
 use discretization
 use geometry_plastic_nonlocal
 use mesh_base
 use DAMASK_interface
 use PETScsys
 use IO
 use debug
 use numerics
 use FEsolving


 implicit none
 private

 integer(pInt), public, protected :: &
   mesh_Nnodes

 integer(pInt), dimension(:), allocatable, private :: &
   microGlobal
 integer(pInt), dimension(:), allocatable, private :: &
   mesh_homogenizationAt

 integer(pInt), dimension(:,:), allocatable, public, protected :: &
   mesh_element                                                                                     !< entryCount and list of elements containing node

 real(pReal), public, protected :: &
   mesh_unitlength                                                                                  !< physical length of one unit in mesh

 real(pReal), dimension(:,:), allocatable, private :: &
   mesh_node                                                                                        !< node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)


 real(pReal), dimension(:,:), allocatable, public, protected :: &
   mesh_ipVolume, &                                                                                 !< volume associated with IP (initially!)
   mesh_node0                                                                                       !< node x,y,z coordinates (initially!)

 real(pReal), dimension(:,:,:), allocatable, public, protected :: &
   mesh_ipArea                                                                                      !< area of interface to neighboring IP (initially!)

 real(pReal), dimension(:,:,:), allocatable, public :: &
   mesh_ipCoordinates                                                                               !< IP x,y,z coordinates (after deformation!)

 real(pReal),dimension(:,:,:,:), allocatable, public, protected :: &
   mesh_ipAreaNormal                                                                                !< area normal of interface to neighboring IP (initially!)

 logical, dimension(3), public, parameter :: mesh_periodicSurface = .true.                          !< flag indicating periodic outer surfaces (used for fluxes)


! grid specific
 integer(pInt), dimension(3), public, protected :: &
   grid                                                                                             !< (global) grid
 integer(pInt), public, protected :: &
   mesh_NcpElemsGlobal, &                                                                           !< total number of CP elements in global mesh
   grid3, &                                                                                         !< (local) grid in 3rd direction
   grid3Offset                                                                                      !< (local) grid offset in 3rd direction
 real(pReal), dimension(3), public, protected :: &
   geomSize
 real(pReal), public, protected :: &
   size3, &                                                                                         !< (local) size in 3rd direction
   size3offset                                                                                      !< (local) size offset in 3rd direction

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
 
  integer(pInt), dimension(3), public :: &
   grid                                                                                             !< (global) grid
 integer(pInt), public :: &
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
 
 call self%tMesh%init('grid',10_pInt,nodes)
 
end subroutine tMesh_grid_init

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

 include 'fftw3-mpi.f03'
 integer(C_INTPTR_T) :: devNull, local_K, local_K_offset
 integer :: ierr, worldsize, j
 integer(pInt), intent(in), optional :: el, ip
 logical :: myDebug

 write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'

 mesh_unitlength = numerics_unitlength                                                              ! set physical extent of a length unit in mesh

 myDebug = (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt)

 call fftw_mpi_init()
 call mesh_spectral_read_grid()


 call MPI_comm_size(PETSC_COMM_WORLD, worldsize, ierr)
 if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_comm_size')
 if(worldsize>grid(3)) call IO_error(894_pInt, ext_msg='number of processes exceeds grid(3)')


 devNull = fftw_mpi_local_size_3d(int(grid(3),C_INTPTR_T), &
                                  int(grid(2),C_INTPTR_T), &
                                  int(grid(1),C_INTPTR_T)/2+1, &
                                  PETSC_COMM_WORLD, &
                                  local_K, &                                                        ! domain grid size along z
                                  local_K_offset)                                                   ! domain grid offset along z
 grid3       = int(local_K,pInt)
 grid3Offset = int(local_K_offset,pInt)
 size3       = geomSize(3)*real(grid3,pReal)      /real(grid(3),pReal)
 size3Offset = geomSize(3)*real(grid3Offset,pReal)/real(grid(3),pReal)

 mesh_NcpElemsGlobal = product(grid)

 mesh_Nnodes  = product(grid(1:2) + 1_pInt)*(grid3 + 1_pInt)

 mesh_node0 = mesh_spectral_build_nodes()
 mesh_node  = mesh_node0
 if (myDebug) write(6,'(a)') ' Built nodes'; flush(6)

 call theMesh%init(mesh_node)
 call theMesh%setNelems(product(grid(1:2))*grid3)
 call mesh_spectral_build_elements()
 mesh_homogenizationAt = mesh_homogenizationAt(product(grid(1:2))*grid3Offset+1: &
                                               product(grid(1:2))*(grid3Offset+grid3))              ! reallocate/shrink in case of MPI
 
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

 if (myDebug) write(6,'(a)') ' Built IP neighborhood'; flush(6)

 if (debug_e < 1 .or. debug_e > theMesh%nElems) &
   call IO_error(602_pInt,ext_msg='element')                                                        ! selected element does not exist
 if (debug_i < 1 .or. debug_i > theMesh%elem%nIPs) &
   call IO_error(602_pInt,ext_msg='IP')                                                             ! selected element does not have requested IP

 FEsolving_execElem = [ 1_pInt,theMesh%nElems ]                                                     ! parallel loop bounds set to comprise all DAMASK elements
 allocate(FEsolving_execIP(2_pInt,theMesh%nElems), source=1_pInt)                                   ! parallel loop bounds set to comprise from first IP...
 forall (j = 1_pInt:theMesh%nElems) FEsolving_execIP(2,j) = theMesh%elem%nIPs                       ! ...up to own IP count for each element


!!!! COMPATIBILITY HACK !!!!
 theMesh%homogenizationAt  = mesh_element(3,:)
 theMesh%microstructureAt  = mesh_element(4,:)
!!!!!!!!!!!!!!!!!!!!!!!!
  call discretization_init(mesh_element(3,:),mesh_element(4,:),&
                           reshape(mesh_ipCoordinates,[3,grid(1)*grid(2)*grid3]),&
                           mesh_node0)

end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Parses geometry file
!> @details important variables have an implicit "save" attribute. Therefore, this function is 
! supposed to be called only once!
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_read_grid()

  character(len=:),            allocatable :: rawData
  character(len=65536)                     :: line
  integer(pInt), allocatable, dimension(:) :: chunkPos
  integer(pInt) :: h =- 1_pInt
  integer(pInt) ::  &
    headerLength = -1_pInt, &                                                                       !< length of header (in lines)
    fileLength, &                                                                                   !< length of the geom file (in characters)
    fileUnit, &
    startPos, endPos, &
    myStat, &
    l, &                                                                                            !< line counter
    c, &                                                                                            !< counter for # microstructures in line
    o, &                                                                                            !< order of "to" packing
    e, &                                                                                            !< "element", i.e. spectral collocation point 
    i, j
    
  grid = -1_pInt
  geomSize = -1.0_pReal

!--------------------------------------------------------------------------------------------------
! read data as stream
  inquire(file = trim(geometryFile), size=fileLength)
  open(newunit=fileUnit, file=trim(geometryFile), access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if(myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=trim(geometryFile))
  allocate(character(len=fileLength)::rawData)
  read(fileUnit) rawData
  close(fileUnit)
  
!--------------------------------------------------------------------------------------------------
! get header length
  endPos = index(rawData,new_line(''))
  if(endPos <= index(rawData,'head')) then
    startPos = len(rawData)
    call IO_error(error_ID=841_pInt, ext_msg='mesh_spectral_read_grid')
  else
    chunkPos = IO_stringPos(rawData(1:endPos))
    if (chunkPos(1) < 2_pInt) call IO_error(error_ID=841_pInt, ext_msg='mesh_spectral_read_grid')
    headerLength = IO_intValue(rawData(1:endPos),chunkPos,1_pInt)
    startPos = endPos + 1_pInt
  endif

!--------------------------------------------------------------------------------------------------
! read and interprete header
  l = 0
  do while (l < headerLength .and. startPos < len(rawData))
    endPos = startPos + index(rawData(startPos:),new_line('')) - 1_pInt
    if (endPos < startPos) endPos = len(rawData)                                                    ! end of file without new line
    line = rawData(startPos:endPos)
    startPos = endPos + 1_pInt
    l = l + 1_pInt

    chunkPos = IO_stringPos(trim(line))
    if (chunkPos(1) < 2) cycle                                                                      ! need at least one keyword value pair
    
    select case ( IO_lc(IO_StringValue(trim(line),chunkPos,1_pInt,.true.)) )
      case ('grid')
        if (chunkPos(1) > 6) then
          do j = 2_pInt,6_pInt,2_pInt
            select case (IO_lc(IO_stringValue(line,chunkPos,j)))
              case('a')
                grid(1) = IO_intValue(line,chunkPos,j+1_pInt)
              case('b')
                grid(2) = IO_intValue(line,chunkPos,j+1_pInt)
              case('c')
                grid(3) = IO_intValue(line,chunkPos,j+1_pInt)
            end select
          enddo
        endif
        
      case ('size')
        if (chunkPos(1) > 6) then
          do j = 2_pInt,6_pInt,2_pInt
            select case (IO_lc(IO_stringValue(line,chunkPos,j)))
              case('x')
                geomSize(1) = IO_floatValue(line,chunkPos,j+1_pInt)
              case('y')
                geomSize(2) = IO_floatValue(line,chunkPos,j+1_pInt)
              case('z')
                geomSize(3) = IO_floatValue(line,chunkPos,j+1_pInt)
            end select
          enddo
        endif
        
      case ('homogenization')
        if (chunkPos(1) > 1) h = IO_intValue(line,chunkPos,2_pInt)
    end select

  enddo

!--------------------------------------------------------------------------------------------------
! sanity checks
  if(h < 1_pInt) &
    call IO_error(error_ID = 842_pInt, ext_msg='homogenization (mesh_spectral_read_grid)')
  if(any(grid < 1_pInt)) &
    call IO_error(error_ID = 842_pInt, ext_msg='grid (mesh_spectral_read_grid)')
  if(any(geomSize < 0.0_pReal)) &
    call IO_error(error_ID = 842_pInt, ext_msg='size (mesh_spectral_read_grid)')

  allocate(microGlobal(product(grid)), source = -1_pInt)
  allocate(mesh_homogenizationAt(product(grid)), source = h)                                        ! too large in case of MPI (shrink later, not very elegant)
     
!--------------------------------------------------------------------------------------------------
! read and interpret content
  e = 1_pInt
  do while (startPos < len(rawData))
    endPos = startPos + index(rawData(startPos:),new_line('')) - 1_pInt
    if (endPos < startPos) endPos = len(rawData)                                                    ! end of file without new line
    line = rawData(startPos:endPos)
    startPos = endPos + 1_pInt
    l = l + 1_pInt
    chunkPos = IO_stringPos(trim(line))
    
    noCompression: if (chunkPos(1) /= 3) then
      c = chunkPos(1)
      microGlobal(e:e+c-1_pInt) =  [(IO_intValue(line,chunkPos,i+1_pInt), i=0_pInt, c-1_pInt)]
    else noCompression
      compression: if (IO_lc(IO_stringValue(line,chunkPos,2))  == 'of') then
        c = IO_intValue(line,chunkPos,1)
        microGlobal(e:e+c-1_pInt) = [(IO_intValue(line,chunkPos,3),i = 1_pInt,IO_intValue(line,chunkPos,1))]
      else if (IO_lc(IO_stringValue(line,chunkPos,2))  == 'to') then compression
        c = abs(IO_intValue(line,chunkPos,3) - IO_intValue(line,chunkPos,1)) + 1_pInt
        o = merge(+1_pInt, -1_pInt, IO_intValue(line,chunkPos,3) > IO_intValue(line,chunkPos,1))
        microGlobal(e:e+c-1_pInt) = [(i, i = IO_intValue(line,chunkPos,1),IO_intValue(line,chunkPos,3),o)]
      else compression
        c = chunkPos(1)
        microGlobal(e:e+c-1_pInt) =  [(IO_intValue(line,chunkPos,i+1_pInt), i=0_pInt, c-1_pInt)]
      endif compression
    endif noCompression

    e = e+c
  end do

  if (e-1 /= product(grid)) call IO_error(error_ID = 843_pInt, el=e)

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
subroutine mesh_spectral_build_elements()

 integer(pInt) :: &
   e, &
   elemOffset

 allocate(mesh_element    (4_pInt+8_pInt,theMesh%nElems), source = 0_pInt)

 elemOffset = product(grid(1:2))*grid3Offset
 do e=1, theMesh%nElems
   mesh_element( 1,e) = -1_pInt                                                                     ! DEPRECATED
   mesh_element( 2,e) = -1_pInt                                                                     ! DEPRECATED
   mesh_element( 3,e) = mesh_homogenizationAt(e)
   mesh_element( 4,e) = microGlobal(e+elemOffset)                                                   ! microstructure
   mesh_element( 5,e) = e + (e-1_pInt)/grid(1) + &
                                     ((e-1_pInt)/(grid(1)*grid(2)))*(grid(1)+1_pInt)                ! base node
   mesh_element( 6,e) = mesh_element(5,e) + 1_pInt
   mesh_element( 7,e) = mesh_element(5,e) + grid(1) + 2_pInt
   mesh_element( 8,e) = mesh_element(5,e) + grid(1) + 1_pInt
   mesh_element( 9,e) = mesh_element(5,e) +(grid(1) + 1_pInt) * (grid(2) + 1_pInt)                  ! second floor base node
   mesh_element(10,e) = mesh_element(9,e) + 1_pInt
   mesh_element(11,e) = mesh_element(9,e) + grid(1) + 2_pInt
   mesh_element(12,e) = mesh_element(9,e) + grid(1) + 1_pInt
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
 integer, dimension(:,:,:,:), allocatable :: &
   ipNeighborhood                                                                                   !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]
 allocate(ipNeighborhood(3,6,1,theMesh%nElems),source=0)

 e = 0_pInt
 do z = 0_pInt,grid3-1_pInt
   do y = 0_pInt,grid(2)-1_pInt
     do x = 0_pInt,grid(1)-1_pInt
       e = e + 1_pInt
        ipNeighborhood(1,1,1,e) = z * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + modulo(x+1_pInt,grid(1)) &
                                      + 1_pInt
         ipNeighborhood(1,2,1,e) = z * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + modulo(x-1_pInt,grid(1)) &
                                      + 1_pInt
         ipNeighborhood(1,3,1,e) = z * grid(1) * grid(2) &
                                      + modulo(y+1_pInt,grid(2)) * grid(1) &
                                      + x &
                                      + 1_pInt
         ipNeighborhood(1,4,1,e) = z * grid(1) * grid(2) &
                                      + modulo(y-1_pInt,grid(2)) * grid(1) &
                                      + x &
                                      + 1_pInt
         ipNeighborhood(1,5,1,e) = modulo(z+1_pInt,grid3) * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + x &
                                      + 1_pInt
         ipNeighborhood(1,6,1,e) = modulo(z-1_pInt,grid3) * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + x &
                                      + 1_pInt
         ipNeighborhood(2,1:6,1,e) = 1_pInt
         ipNeighborhood(3,1,1,e) = 2_pInt
         ipNeighborhood(3,2,1,e) = 1_pInt
         ipNeighborhood(3,3,1,e) = 4_pInt
         ipNeighborhood(3,4,1,e) = 3_pInt
         ipNeighborhood(3,5,1,e) = 6_pInt
         ipNeighborhood(3,6,1,e) = 5_pInt
     enddo
   enddo
 enddo
 
 call geometry_plastic_nonlocal_set_IPneighborhood(ipNeighborhood)

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

 integer(pInt) :: &
   i,j,k,n
 integer(pInt),           dimension(3), parameter :: &
   diag = 1_pInt
 integer(pInt),           dimension(3) :: &
   shift = 0_pInt, &
   lookup = 0_pInt, &
   me = 0_pInt, &
   iRes = 0_pInt
 integer(pInt),           dimension(3,8) :: &
   neighbor = reshape([ &
                       0_pInt, 0_pInt, 0_pInt, &
                       1_pInt, 0_pInt, 0_pInt, &
                       1_pInt, 1_pInt, 0_pInt, &
                       0_pInt, 1_pInt, 0_pInt, &
                       0_pInt, 0_pInt, 1_pInt, &
                       1_pInt, 0_pInt, 1_pInt, &
                       1_pInt, 1_pInt, 1_pInt, &
                       0_pInt, 1_pInt, 1_pInt  ], [3,8])

!--------------------------------------------------------------------------------------------------
! initializing variables
 iRes =  [size(centres,2),size(centres,3),size(centres,4)]
 nodes = 0.0_pReal
 wrappedCentres = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! report
 if (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt) then
   write(6,'(a)')          ' Meshing cubes around centroids'
   write(6,'(a,3(e12.5))') ' Dimension: ', gDim
   write(6,'(a,3(i5))')    ' Resolution:', iRes
 endif

!--------------------------------------------------------------------------------------------------
! building wrappedCentres = centroids + ghosts
 wrappedCentres(1:3,2_pInt:iRes(1)+1_pInt,2_pInt:iRes(2)+1_pInt,2_pInt:iRes(3)+1_pInt) = centres
 do k = 0_pInt,iRes(3)+1_pInt
   do j = 0_pInt,iRes(2)+1_pInt
     do i = 0_pInt,iRes(1)+1_pInt
       if (k==0_pInt .or. k==iRes(3)+1_pInt .or. &                                                  ! z skin
           j==0_pInt .or. j==iRes(2)+1_pInt .or. &                                                  ! y skin
           i==0_pInt .or. i==iRes(1)+1_pInt      ) then                                             ! x skin
         me = [i,j,k]                                                                               ! me on skin
         shift = sign(abs(iRes+diag-2_pInt*me)/(iRes+diag),iRes+diag-2_pInt*me)
         lookup = me-diag+shift*iRes
         wrappedCentres(1:3,i+1_pInt,        j+1_pInt,        k+1_pInt) = &
                centres(1:3,lookup(1)+1_pInt,lookup(2)+1_pInt,lookup(3)+1_pInt) &
                - matmul(Favg, real(shift,pReal)*gDim)
       endif
 enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! averaging
 do k = 0_pInt,iRes(3); do j = 0_pInt,iRes(2); do i = 0_pInt,iRes(1)
   do n = 1_pInt,8_pInt
    nodes(1:3,i+1_pInt,j+1_pInt,k+1_pInt) = &
    nodes(1:3,i+1_pInt,j+1_pInt,k+1_pInt) + wrappedCentres(1:3,i+1_pInt+neighbor(1,n), &
                                                               j+1_pInt+neighbor(2,n), &
                                                               k+1_pInt+neighbor(3,n) )
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
