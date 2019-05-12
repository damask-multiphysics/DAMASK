!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc, Abaqus and the spectral solver
!--------------------------------------------------------------------------------------------------
module mesh
 use, intrinsic :: iso_c_binding
 use prec, only: pReal, pInt
 use mesh_base

 implicit none
 private
 integer(pInt), public, protected :: &
   mesh_Nnodes, &                                                                                   !< total number of nodes in mesh
   mesh_Ncellnodes, &                                                                               !< total number of cell nodes in mesh (including duplicates)
   mesh_Ncells, &                                                                                   !< total number of cells in mesh
   mesh_maxNipNeighbors, &                                                                          !< max number of IP neighbors in any CP element
   mesh_maxNsharedElems                                                                             !< max number of CP elements sharing a node


 integer(pInt), dimension(:), allocatable, private :: &
   microGlobal
 integer(pInt), dimension(:), allocatable, private :: &
   mesh_homogenizationAt

 integer(pInt), dimension(:,:), allocatable, public, protected :: &
   mesh_element                                                                                     !< entryCount and list of elements containing node

 integer(pInt), dimension(:,:,:,:), allocatable, public, protected :: &
   mesh_ipNeighborhood                                                                              !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]

 real(pReal), public, protected :: &
   mesh_unitlength                                                                                  !< physical length of one unit in mesh

 real(pReal), dimension(:,:), allocatable, public :: &
   mesh_node, &                                                                                     !< node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)
   mesh_cellnode                                                                                    !< cell node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)

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

integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_cellnodeParent                                                                              !< cellnode's parent element ID, cellnode's intra-element ID

 integer(pInt),dimension(:,:,:), allocatable, private :: &
   mesh_cell                                                                                        !< cell connectivity for each element,ip/cell

 integer(pInt), dimension(:,:,:), allocatable, private :: &
   FE_cellface                                                                                      !< list of intra-cell cell node IDs that constitute the cell faces of a specific type of cell


! These definitions should actually reside in the FE-solver specific part (different for MARC/ABAQUS)
! Hence, I suggest to prefix with "FE_"

 integer(pInt), parameter, private :: &
   FE_Ngeomtypes = 10_pInt, &
   FE_Ncelltypes = 4_pInt, &
   FE_maxNmatchingNodesPerFace = 4_pInt, &
   FE_maxNfaces = 6_pInt, &
   FE_maxNcellnodesPerCell = 8_pInt, &
   FE_maxNcellfaces = 6_pInt, &
   FE_maxNcellnodesPerCellface = 4_pInt



 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCell = &             !< number of cell nodes in a specific cell type
 int([ &
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      8  & ! (3D 8node)
  ],pInt)

 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCellface = &        !< number of cell nodes per cell face in a specific cell type
 int([&
      2, & ! (2D 3node)
      2, & ! (2D 4node)
      3, & ! (3D 4node)
      4  & ! (3D 8node)
  ],pInt)


 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NipNeighbors = &                  !< number of ip neighbors / cell faces in a specific cell type
 int([&
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      6  & ! (3D 8node)
  ],pInt)


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
   mesh_init, &
   mesh_cellCenterCoordinates

 private :: &
   mesh_build_cellconnectivity, &
   mesh_build_ipAreas, &
   mesh_build_FEdata, &
   mesh_spectral_build_nodes, &
   mesh_spectral_build_elements, &
   mesh_spectral_build_ipNeighborhood, &
   mesh_build_cellnodes, &
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
 
 implicit none
 class(tMesh_grid) :: self
 real(pReal), dimension(:,:), intent(in) :: nodes
 
 call self%tMesh%init('grid',10_pInt,nodes)
 
end subroutine tMesh_grid_init

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

#include <petsc/finclude/petscsys.h>
 use PETScsys

 use DAMASK_interface
 use IO, only: &
   IO_error
 use debug, only: &
   debug_e, &
   debug_i, &
   debug_level, &
   debug_mesh, &
   debug_levelBasic
 use numerics, only: &
   numerics_unitlength
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP

 implicit none
 include 'fftw3-mpi.f03'
 integer(C_INTPTR_T) :: devNull, local_K, local_K_offset
 integer :: ierr, worldsize, i
 integer(pInt), intent(in), optional :: el, ip
 integer(pInt) :: j
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

 call mesh_spectral_build_nodes()
 if (myDebug) write(6,'(a)') ' Built nodes'; flush(6)

 call theMesh%init(mesh_node)
 call theMesh%setNelems(product(grid(1:2))*grid3)
 mesh_homogenizationAt = mesh_homogenizationAt(product(grid(1:2))*grid3)                            ! reallocate/shrink in case of MPI
 mesh_maxNipNeighbors = theMesh%elem%nIPneighbors
 
 call mesh_spectral_build_elements()
 if (myDebug) write(6,'(a)') ' Built elements'; flush(6)
 
 
 call mesh_build_FEdata                                                                             ! get properties of the different types of elements
 call mesh_build_cellconnectivity
 if (myDebug) write(6,'(a)') ' Built cell connectivity'; flush(6)
 mesh_cellnode = mesh_build_cellnodes(mesh_node,mesh_Ncellnodes)
 if (myDebug) write(6,'(a)') ' Built cell nodes'; flush(6)
 call mesh_build_ipCoordinates
 if (myDebug) write(6,'(a)') ' Built IP coordinates'; flush(6)
 allocate(mesh_ipVolume(1,theMesh%nElems),source=product([geomSize(1:2),size3]/real([grid(1:2),grid3])))
 if (myDebug) write(6,'(a)') ' Built IP volumes'; flush(6)
 call mesh_build_ipAreas
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
! for a homogeneous mesh, all elements have the same number of IPs and and cell nodes.
! hence, xxPerElem instead of maxXX
! better name
 theMesh%homogenizationAt  = mesh_element(3,:)
 theMesh%microstructureAt  = mesh_element(4,:)
!!!!!!!!!!!!!!!!!!!!!!!!
 deallocate(mesh_cell)
end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Parses geometry file
!> @details important variables have an implicit "save" attribute. Therefore, this function is 
! supposed to be called only once!
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_read_grid()
 use IO, only: &
   IO_stringPos, &
   IO_lc, &
   IO_stringValue, &
   IO_intValue, &
   IO_floatValue, &
   IO_error
 use DAMASK_interface, only: &
   geometryFile

  implicit none
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
! read and interprete content
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


!--------------------------------------------------------------------------------------------------
!> @brief Store x,y,z coordinates of all nodes in mesh.
!! Allocates global arrays 'mesh_node0' and 'mesh_node'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_nodes()

 implicit none
 integer(pInt) :: n

 allocate (mesh_node0 (3,mesh_Nnodes), source = 0.0_pReal)

 forall (n = 0_pInt:mesh_Nnodes-1_pInt)
   mesh_node0(1,n+1_pInt) = mesh_unitlength * &
           geomSize(1)*real(mod(n,(grid(1)+1_pInt) ),pReal) &
                                                  / real(grid(1),pReal)
   mesh_node0(2,n+1_pInt) = mesh_unitlength * &
           geomSize(2)*real(mod(n/(grid(1)+1_pInt),(grid(2)+1_pInt)),pReal) &
                                                  / real(grid(2),pReal)
   mesh_node0(3,n+1_pInt) = mesh_unitlength * &
           size3*real(mod(n/(grid(1)+1_pInt)/(grid(2)+1_pInt),(grid3+1_pInt)),pReal) &
                                                  / real(grid3,pReal) + &
           size3offset
 end forall

 mesh_node = mesh_node0

end subroutine mesh_spectral_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, material, texture, and node list per element.
!! Allocates global array 'mesh_element'
!> @todo does the IO_error makes sense?
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_elements()
 use IO, only: &
   IO_error
 implicit none
 integer(pInt) :: &
   e, &
   elemOffset


 allocate(mesh_element    (4_pInt+8_pInt,theMesh%nElems), source = 0_pInt)

 elemOffset = product(grid(1:2))*grid3Offset
 e = 0_pInt
 do while (e < theMesh%nElems)                                                                      ! fill expected number of elements, stop at end of data
   e = e+1_pInt                                                                                     ! valid element entry
   mesh_element( 1,e) = -1_pInt                                                                     ! DEPRECATED
   mesh_element( 2,e) = 10_pInt
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

 if (e /= theMesh%nElems) call IO_error(880_pInt,e)

end subroutine mesh_spectral_build_elements


!--------------------------------------------------------------------------------------------------
!> @brief build neighborhood relations for spectral
!> @details assign globals: mesh_ipNeighborhood
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_ipNeighborhood

 implicit none
 integer(pInt) :: &
  x,y,z, &
  e
 allocate(mesh_ipNeighborhood(3,6,1,theMesh%nElems),source=0_pInt)

 e = 0_pInt
 do z = 0_pInt,grid3-1_pInt
   do y = 0_pInt,grid(2)-1_pInt
     do x = 0_pInt,grid(1)-1_pInt
       e = e + 1_pInt
         mesh_ipNeighborhood(1,1,1,e) = z * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + modulo(x+1_pInt,grid(1)) &
                                      + 1_pInt
         mesh_ipNeighborhood(1,2,1,e) = z * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + modulo(x-1_pInt,grid(1)) &
                                      + 1_pInt
         mesh_ipNeighborhood(1,3,1,e) = z * grid(1) * grid(2) &
                                      + modulo(y+1_pInt,grid(2)) * grid(1) &
                                      + x &
                                      + 1_pInt
         mesh_ipNeighborhood(1,4,1,e) = z * grid(1) * grid(2) &
                                      + modulo(y-1_pInt,grid(2)) * grid(1) &
                                      + x &
                                      + 1_pInt
         mesh_ipNeighborhood(1,5,1,e) = modulo(z+1_pInt,grid3) * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + x &
                                      + 1_pInt
         mesh_ipNeighborhood(1,6,1,e) = modulo(z-1_pInt,grid3) * grid(1) * grid(2) &
                                      + y * grid(1) &
                                      + x &
                                      + 1_pInt
         mesh_ipNeighborhood(2,1:6,1,e) = 1_pInt
         mesh_ipNeighborhood(3,1,1,e) = 2_pInt
         mesh_ipNeighborhood(3,2,1,e) = 1_pInt
         mesh_ipNeighborhood(3,3,1,e) = 4_pInt
         mesh_ipNeighborhood(3,4,1,e) = 3_pInt
         mesh_ipNeighborhood(3,5,1,e) = 6_pInt
         mesh_ipNeighborhood(3,6,1,e) = 5_pInt
     enddo
   enddo
 enddo

end subroutine mesh_spectral_build_ipNeighborhood


!--------------------------------------------------------------------------------------------------
!> @brief builds mesh of (distorted) cubes for given coordinates (= center of the cubes)
!--------------------------------------------------------------------------------------------------
function mesh_nodesAroundCentres(gDim,Favg,centres) result(nodes)
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic

 implicit none
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


!#################################################################################################################
!#################################################################################################################
!#################################################################################################################
! The following routines are not solver specific and should be included in mesh_base (most likely in modified form)
!#################################################################################################################
!#################################################################################################################
!#################################################################################################################



!--------------------------------------------------------------------------------------------------
!> @brief Split CP elements into cells.
!> @details Build a mapping between cells and the corresponding cell nodes ('mesh_cell').
!> Cell nodes that are also matching nodes are unique in the list of cell nodes,
!> all others (currently) might be stored more than once.
!> Also allocates the 'mesh_node' array.
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_cellconnectivity

 implicit none
 integer(pInt), dimension(:), allocatable :: &
   matchingNode2cellnode
 integer(pInt), dimension(:,:), allocatable :: &
   cellnodeParent
 integer(pInt), dimension(theMesh%elem%Ncellnodes) :: &
   localCellnode2globalCellnode
 integer(pInt) :: &
   e,n,i, &
   matchingNodeID, &
   localCellnodeID
   
 integer(pInt), dimension(FE_Ngeomtypes), parameter :: FE_NmatchingNodes = &               !< number of nodes that are needed for face matching in a specific type of element geometry
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      3, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      4, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      8, & ! element   7 (3D 8node 8ip)
      8  & ! element  21 (3D 20node 27ip)
  ],pInt)

 allocate(mesh_cell(FE_maxNcellnodesPerCell,theMesh%elem%nIPs,theMesh%nElems), source=0_pInt)
 allocate(matchingNode2cellnode(theMesh%nNodes),                            source=0_pInt)
 allocate(cellnodeParent(2_pInt,theMesh%elem%Ncellnodes*theMesh%nElems),       source=0_pInt)
 
 mesh_Ncells = theMesh%nElems*theMesh%elem%nIPs
!--------------------------------------------------------------------------------------------------
! Count cell nodes (including duplicates) and generate cell connectivity list
 mesh_Ncellnodes = 0_pInt

 do e = 1_pInt,theMesh%nElems
   localCellnode2globalCellnode = 0_pInt
   do i = 1_pInt,theMesh%elem%nIPs
     do n = 1_pInt,theMesh%elem%NcellnodesPerCell
       localCellnodeID = theMesh%elem%cell(n,i)
       if (localCellnodeID <= FE_NmatchingNodes(theMesh%elem%geomType)) then                                            ! this cell node is a matching node
         matchingNodeID = mesh_element(4_pInt+localCellnodeID,e)
         if (matchingNode2cellnode(matchingNodeID) == 0_pInt) then                                  ! if this matching node does not yet exist in the glbal cell node list ...
           mesh_Ncellnodes = mesh_Ncellnodes + 1_pInt                                               ! ... count it as cell node ...
           matchingNode2cellnode(matchingNodeID) = mesh_Ncellnodes                                  ! ... and remember its global ID
           cellnodeParent(1_pInt,mesh_Ncellnodes) = e                                               ! ... and where it belongs to
           cellnodeParent(2_pInt,mesh_Ncellnodes) = localCellnodeID
         endif
         mesh_cell(n,i,e) = matchingNode2cellnode(matchingNodeID)
       else                                                                                         ! this cell node is no matching node
         if (localCellnode2globalCellnode(localCellnodeID) == 0_pInt) then                          ! if this local cell node does not yet exist in the  global cell node list ...
           mesh_Ncellnodes = mesh_Ncellnodes + 1_pInt                                               ! ... count it as cell node ...
           localCellnode2globalCellnode(localCellnodeID) = mesh_Ncellnodes                          ! ... and remember its global ID ...
           cellnodeParent(1_pInt,mesh_Ncellnodes) = e                                               ! ... and it belongs to
           cellnodeParent(2_pInt,mesh_Ncellnodes) = localCellnodeID
         endif
         mesh_cell(n,i,e) = localCellnode2globalCellnode(localCellnodeID)
       endif
     enddo
   enddo
 enddo

 allocate(mesh_cellnodeParent(2_pInt,mesh_Ncellnodes))
 allocate(mesh_cellnode(3_pInt,mesh_Ncellnodes))
 
 forall(n = 1_pInt:mesh_Ncellnodes)
   mesh_cellnodeParent(1,n) = cellnodeParent(1,n)
   mesh_cellnodeParent(2,n) = cellnodeParent(2,n)
 endforall

end subroutine mesh_build_cellconnectivity


!--------------------------------------------------------------------------------------------------
!> @brief Calculate position of cellnodes from the given position of nodes
!> Build list of cellnodes' coordinates.
!> Cellnode coordinates are calculated from a weighted sum of node coordinates.
!--------------------------------------------------------------------------------------------------
function mesh_build_cellnodes(nodes,Ncellnodes)

 implicit none
 integer(pInt),                         intent(in) :: Ncellnodes                                    !< requested number of cellnodes
 real(pReal), dimension(3,mesh_Nnodes), intent(in) :: nodes
 real(pReal), dimension(3,Ncellnodes) :: mesh_build_cellnodes

 integer(pInt) :: &
   e,n,m, &
   localCellnodeID
 real(pReal), dimension(3) :: &
   myCoords

 mesh_build_cellnodes = 0.0_pReal
!$OMP PARALLEL DO PRIVATE(e,localCellnodeID,myCoords)
 do n = 1_pInt,Ncellnodes                                                                           ! loop over cell nodes
   e = mesh_cellnodeParent(1,n)
   localCellnodeID = mesh_cellnodeParent(2,n)
   myCoords = 0.0_pReal
   do m = 1_pInt,theMesh%elem%nNodes
     myCoords = myCoords + nodes(1:3,mesh_element(4_pInt+m,e)) &
                         * theMesh%elem%cellNodeParentNodeWeights(m,localCellnodeID)
   enddo
   mesh_build_cellnodes(1:3,n) = myCoords / sum(theMesh%elem%cellNodeParentNodeWeights(:,localCellnodeID))
 enddo
!$OMP END PARALLEL DO

end function mesh_build_cellnodes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP Coordinates. Allocates global array 'mesh_ipCoordinates'
! Called by all solvers in mesh_init in order to initialize the ip coordinates.
! Later on the current ip coordinates are directly prvided by the spectral solver and by Abaqus,
! so no need to use this subroutine anymore; Marc however only provides nodal displacements,
! so in this case the ip coordinates are always calculated on the basis of this subroutine.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR THE MOMENT THIS SUBROUTINE ACTUALLY CALCULATES THE CELL CENTER AND NOT THE IP COORDINATES,
! AS THE IP IS NOT (ALWAYS) LOCATED IN THE CENTER OF THE IP VOLUME.
! HAS TO BE CHANGED IN A LATER VERSION.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipCoordinates

 implicit none
 integer(pInt) :: e,c,i,n
 real(pReal), dimension(3) :: myCoords

 allocate(mesh_ipCoordinates(3,theMesh%elem%nIPs,theMesh%nElems),source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(c,myCoords)
 do e = 1_pInt,theMesh%nElems                                                                        ! loop over cpElems
   c = theMesh%elem%cellType
   do i = 1_pInt,theMesh%elem%nIPs                                                                         ! loop over ips=cells in this element
     myCoords = 0.0_pReal
     do n = 1_pInt,FE_NcellnodesPerCell(c)                                                          ! loop over cell nodes in this cell
       myCoords = myCoords + mesh_cellnode(1:3,mesh_cell(n,i,e))
     enddo
     mesh_ipCoordinates(1:3,i,e) = myCoords / real(FE_NcellnodesPerCell(c),pReal)
   enddo
 enddo
 !$OMP END PARALLEL DO

end subroutine mesh_build_ipCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief Calculates cell center coordinates.
!--------------------------------------------------------------------------------------------------
pure function mesh_cellCenterCoordinates(ip,el)

 implicit none
 integer(pInt), intent(in) :: el, &                                                                  !< element number
                              ip                                                                     !< integration point number
 real(pReal), dimension(3) :: mesh_cellCenterCoordinates                                             !< x,y,z coordinates of the cell center of the requested IP cell
 integer(pInt) :: c,n

 c = theMesh%elem%cellType
 mesh_cellCenterCoordinates = 0.0_pReal
 do n = 1_pInt,FE_NcellnodesPerCell(c)                                                               ! loop over cell nodes in this cell
   mesh_cellCenterCoordinates = mesh_cellCenterCoordinates + mesh_cellnode(1:3,mesh_cell(n,ip,el))
 enddo
 mesh_cellCenterCoordinates = mesh_cellCenterCoordinates / real(FE_NcellnodesPerCell(c),pReal)

end function mesh_cellCenterCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipAreas
 use math, only: &
   math_crossproduct

 implicit none
 integer(pInt) :: e,t,g,c,i,f,n,m
 real(pReal), dimension (3,FE_maxNcellnodesPerCellface) :: nodePos, normals
 real(pReal), dimension(3) :: normal

 allocate(mesh_ipArea(theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems), source=0.0_pReal)
 allocate(mesh_ipAreaNormal(3_pInt,theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems), source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,nodePos,normal,normals)
   do e = 1_pInt,theMesh%nElems                                                                      ! loop over cpElems
     c = theMesh%elem%cellType
     select case (c)

       case (1_pInt,2_pInt)                                                                         ! 2D 3 or 4 node
         do i = 1_pInt,theMesh%elem%nIPs                                                                   ! loop over ips=cells in this element
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal(1) =   nodePos(2,2) - nodePos(2,1)                                              ! x_normal =  y_connectingVector
             normal(2) = -(nodePos(1,2) - nodePos(1,1))                                             ! y_normal = -x_connectingVector
             normal(3) = 0.0_pReal
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (3_pInt)                                                                                ! 3D 4node
         do i = 1_pInt,theMesh%elem%nIPs                                                                   ! loop over ips=cells in this element
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal = math_crossproduct(nodePos(1:3,2) - nodePos(1:3,1), &
                                         nodePos(1:3,3) - nodePos(1:3,1))
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (4_pInt)                                                                                ! 3D 8node
         ! for this cell type we get the normal of the quadrilateral face as an average of
         ! four normals of triangular subfaces; since the face consists only of two triangles,
         ! the sum has to be divided by two; this whole prcedure tries to compensate for
         ! probable non-planar cell surfaces
         m = FE_NcellnodesPerCellface(c)
         do i = 1_pInt,theMesh%elem%nIPs                                                                   ! loop over ips=cells in this element
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               normals(1:3,n) = 0.5_pReal &
                              * math_crossproduct(nodePos(1:3,1+mod(n  ,m)) - nodePos(1:3,n), &
                                                   nodePos(1:3,1+mod(n+1,m)) - nodePos(1:3,n))
             normal = 0.5_pReal * sum(normals,2)
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)
           enddo
         enddo

     end select
   enddo
 !$OMP END PARALLEL DO

end subroutine mesh_build_ipAreas


!--------------------------------------------------------------------------------------------------
!> @brief get properties of different types of finite elements
!> @details assign globals: FE_nodesAtIP, FE_ipNeighbor, FE_subNodeOnIPFace
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_FEdata

 implicit none
 integer(pInt) :: me
 allocate(FE_cellface(FE_maxNcellnodesPerCellface,FE_maxNcellfaces,FE_Ncelltypes),  source=0_pInt)


 ! *** FE_cellface ***
 me = 0_pInt

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 2D 3node, VTK_TRIANGLE (5)
    reshape(int([&
    2,3,  &
    3,1,  &
    1,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 2D 4node, VTK_QUAD (9)
    reshape(int([&
    2,3,  &
    4,1,  &
    3,4,  &
    1,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 3D 4node, VTK_TETRA (10)
    reshape(int([&
    1,3,2,  &
    1,2,4,  &
    2,3,4,  &
    1,4,3   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 3D 8node, VTK_HEXAHEDRON (12)
    reshape(int([&
    2,3,7,6,  &
    4,1,5,8,  &
    3,4,8,7,  &
    1,2,6,5,  &
    5,6,7,8,  &
    1,4,3,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])
  

end subroutine mesh_build_FEdata


end module mesh
