! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------------------------------
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!! Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!! Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!! Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!! Krishna Komerla, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc, Abaqus and the spectral solver 
!--------------------------------------------------------------------------------------------------

module mesh     
 use, intrinsic :: iso_c_binding
 use prec, only: pReal, pInt

 implicit none
 private
 
 integer(pInt), public, protected :: &
   mesh_NcpElems, &                                                                                 !< total number of CP elements in mesh
   mesh_NelemSets, &
   mesh_maxNelemInSet, &
   mesh_Nmaterials, &
   mesh_Nnodes, &                                                                                   !< total number of nodes in mesh
   mesh_maxNnodes, &                                                                                !< max number of nodes in any CP element
   mesh_maxNips, &                                                                                  !< max number of IPs in any CP element
   mesh_maxNipNeighbors, &                                                                          !< max number of IP neighbors in any CP element
   mesh_maxNsharedElems, &                                                                          !< max number of CP elements sharing a node
   mesh_maxNsubNodes, &
   mesh_Nelems                                                                                      !< total number of elements in mesh

 integer(pInt), dimension(:,:), allocatable, public, protected :: &
   mesh_element, &                                                                                  !< FEid, type(internal representation), material, texture, node indices
   mesh_sharedElem, &                                                                               !< entryCount and list of elements containing node
   mesh_nodeTwins                                                                                   !< node twins are surface nodes that lie exactly on opposite sides of the mesh (surfaces nodes with equal coordinate values in two dimensions)
 
 integer(pInt), dimension(:,:,:,:), allocatable, public, protected :: &
   mesh_ipNeighborhood                                                                              !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]

 real(pReal), dimension(:,:), allocatable, public :: &
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
    
 logical, dimension(3), public, protected :: mesh_periodicSurface                                   !< flag indicating periodic outer surfaces (used for fluxes)

#ifdef Marc
 integer(pInt), private :: &  
   hypoelasticTableStyle, &                                                                         !< Table style (Marc only)
   initialcondTableStyle                                                                            !< Table style (Marc only)
#endif
 
 integer(pInt), dimension(2), private :: &
   mesh_maxValStateVar = 0_pInt
             
 character(len=64), dimension(:), allocatable, private :: &
   mesh_nameElemSet, &                                                                              !< names of elementSet
   mesh_nameMaterial, &                                                                             !< names of material in solid section
   mesh_mapMaterial                                                                                 !< name of elementSet for material
     
 integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
     
 integer(pInt), dimension(:,:), allocatable, target, private :: &
   mesh_mapFEtoCPelem, &                                                                            !< [sorted FEid, corresponding CPid]
   mesh_mapFEtoCPnode                                                                               !< [sorted FEid, corresponding CPid]
   
 real(pReal),dimension(:,:,:), allocatable, private :: &
   mesh_subNodeCoord                                                                                !< coordinates of subnodes per element
 
 logical, private :: noPart                                                                         !< for cases where the ABAQUS input file does not use part/assembly information

#ifdef Spectral
 include 'fftw3.f03'
 real(pReal),   dimension(3), public, protected :: &
   geomdim, &                                                                                       !< physical dimension of volume element per direction
   scaledDim                                                                                        !< scaled dimension of volume element, depending on selected divergence calculation
 integer(pInt), dimension(3), public, protected :: &
   res                                                                                              !< resolution, e.g. number of Fourier points in each direction
 real(pReal),   public, protected  :: &
   wgt
 integer(pInt), public, protected  :: &
   res1_red, &
   homog
#endif

! These definitions should actually reside in the FE-solver specific part (different for MARC/ABAQUS)
! Hence, I suggest to prefix with "FE_"

 integer(pInt), parameter, public :: &
   FE_Nelemtypes       = 13_pInt, &
   FE_Ngeomtypes       = 10_pInt, &
   FE_maxNnodes        = 8_pInt, &
   FE_maxNsubNodes     = 56_pInt, &
   FE_maxNips          = 27_pInt, &
   FE_maxNipNeighbors  = 6_pInt, &
   FE_maxmaxNnodesAtIP = 8_pInt, &                                                                  !< max number of (equivalent) nodes attached to an IP
   FE_maxNipFaceNodes  = 4_pInt
                      
 integer(pInt), dimension(FE_Nelemtypes), parameter, public :: FE_geomtype = &                      !< geometry type of particular element type
 int([ &
      1, & ! element   6 (2D 3node 1ip)
      2, & ! element 125 (2D 6node 3ip)
      3, & ! element  11 (2D 4node 4ip)
      4, & ! element  27 (2D 8node 9ip)
      3, & ! element  54 (2D 8node 4ip)
      5, & ! element 134 (3D 4node 1ip)
      6, & ! element 157 (3D 5node 4ip)
      6, & ! element 127 (3D 10node 4ip)
      7, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      9, & ! element   7 (3D 8node 8ip)
      9, & ! element  57 (3D 20node 8ip)
     10  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Nelemtypes), parameter, private :: FE_NoriginalNodes = &               !< nodes in a specific type of element (how it is originally defined by marc)
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      6, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      8, & ! element  27 (2D 8node 9ip)
      8, & ! element  54 (2D 8node 4ip)
      4, & ! element 134 (3D 4node 1ip)
      5, & ! element 157 (3D 5node 4ip)
     10, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      8, & ! element   7 (3D 8node 8ip)
     20, & ! element  57 (3D 20node 8ip)
     20  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_dimension = &                        !< dimension of element type
 int([ &
      2, & ! element   6 (2D 3node 1ip)
      2, & ! element 125 (2D 6node 3ip)
      2, & ! element  11 (2D 4node 4ip)
      2, & ! element  27 (2D 8node 9ip)
      3, & ! element 134 (3D 4node 1ip)
      3, & ! element 127 (3D 10node 4ip)
      3, & ! element 136 (3D 6node 6ip)
      3, & ! element 117 (3D 8node 1ip)
      3, & ! element   7 (3D 8node 8ip)
      3  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_Nnodes = &                        !< nodes in a specific type of element (how we use it) 
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

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_Nips = &                          !< IPs in a specific type of element
 int([ &
      1, & ! element   6 (2D 3node 1ip)
      3, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      9, & ! element  27 (2D 8node 9ip)
      1, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      1, & ! element 117 (3D 8node 1ip)
      8, & ! element   7 (3D 8node 8ip)
     27  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_NipNeighbors = &                  !< IP neighbors in a specific type of element
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      4, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      4, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
      6, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      6, & ! element 117 (3D 8node 1ip)
      6, & ! element   7 (3D 8node 8ip)
      6  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, private  :: FE_NsubNodes = &                   !< subnodes required to fully define all IP volumes
 int([ &
      0, & ! element   6 (2D 3node 1ip)
      4, & ! element 125 (2D 6node 3ip)
      5, & ! element  11 (2D 4node 4ip)
     12, & ! element  27 (2D 8node 9ip)
      0, & ! element 134 (3D 4node 1ip)
     11, & ! element 127 (3D 10node 4ip)
     15, & ! element 136 (3D 6node 6ip)
      0, & ! element 117 (3D 8node 1ip)
     19, & ! element   7 (3D 8node 8ip)
     56  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_maxNipNeighbors,FE_Ngeomtypes), parameter, private :: FE_NfaceNodes = &!< nodes per face in a specific type of element
 reshape(int([ &
  2,2,2,0,0,0, & ! element   6 (2D 6node 1ip)
  2,2,2,0,0,0, & ! element 125 (2D 6node 3ip)
  2,2,2,2,0,0, & ! element  11 (2D 4node 4ip)
  2,2,2,2,0,0, & ! element  27 (2D 8node 9ip)
  3,3,3,3,0,0, & ! element 134 (3D 4node 1ip)
  3,3,3,3,0,0, & ! element 127 (3D 10node 4ip)
  3,4,4,4,3,0, & ! element 136 (3D 6node 6ip)
  4,4,4,4,4,4, & ! element 117 (3D 8node 1ip)
  4,4,4,4,4,4, & ! element   7 (3D 8node 8ip)
  4,4,4,4,4,4  & ! element  21 (3D 20node 27ip)
  ],pInt),[FE_maxNipNeighbors,FE_Ngeomtypes])   

 integer(pInt), dimension(FE_Ngeomtypes), parameter, private :: FE_maxNnodesAtIP = &                !< map IP index to two node indices in a specific type of element
 int([ &
      3, & ! element   6 (2D 6node 1ip)
      1, & ! element 125 (2D 6node 3ip)
      1, & ! element  11 (2D 4node 4ip)
      2, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
      1, & ! element 127 (3D 10node 4ip)
      1, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      1, & ! element   7 (3D 8node 8ip)
      4  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, private :: FE_NipFaceNodes = &                !< number of subnodes that constitute an ip face of a specific type of element
 int([ &
      2, & ! element   6 (2D 6node 1ip)
      2, & ! element 125 (2D 6node 3ip)
      2, & ! element  11 (2D 4node 4ip)
      2, & ! element  27 (2D 8node 9ip)
      3, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      4, & ! element 136 (3D 6node 6ip)
      4, & ! element 117 (3D 8node 1ip)
      4, & ! element   7 (3D 8node 8ip)
      4  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_maxNipFaceNodes,FE_maxNipNeighbors,FE_Ngeomtypes), parameter, private :: &
                                                                           FE_nodeOnFace = &        !< List of node indices on each face of a specific type of element
 reshape(int([&
  1,2,0,0 , & ! element   6 (2D 3node 1ip)
  2,3,0,0 , &
  3,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element 125 (2D 6node 3ip)
  2,3,0,0 , &
  3,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element  11 (2D 4node 4ip)
  2,3,0,0 , &
  3,4,0,0 , &
  4,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element  27 (2D 8node 9ip)
  2,3,0,0 , &
  3,4,0,0 , &
  4,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 134 (3D 4node 1ip)
  1,4,2,0 , &
  2,3,4,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 127 (3D 10node 4ip)
  1,4,2,0 , &
  2,4,3,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 136 (3D 6node 6ip)
  1,4,5,2 , &
  2,5,6,3 , &
  1,3,6,4 , &
  4,6,5,0 , &
  0,0,0,0 , &
  1,2,3,4 , & ! element 117 (3D 8node 1ip)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,4 , & ! element   7 (3D 8node 8ip)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,4 , & ! element  21 (3D 20node 27ip)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5   &
  ],pInt),[FE_maxNipFaceNodes,FE_maxNipNeighbors,FE_Ngeomtypes])
 
 integer(pInt), dimension(:,:,:), allocatable, private :: &
   FE_nodesAtIP, &                                                                                  !< map IP index to two node indices in a specific type of element
   FE_ipNeighbor, &                                                                                 !< +x,-x,+y,-y,+z,-z list of intra-element IPs and(negative) neighbor faces per own IP in a specific type of element
   FE_subNodeParent
  
 integer(pInt), dimension(:,:,:,:), allocatable, private :: &
   FE_subNodeOnIPFace
       
 public  :: mesh_init, &
            mesh_FEasCP, &
            mesh_build_subNodeCoords, &
            mesh_build_ipVolumes, &
            mesh_build_ipCoordinates, &
            mesh_cellCenterCoordinates
#ifdef Spectral
 public  :: mesh_regrid, &
            mesh_nodesAroundCentres, &
            mesh_deformedCoordsFFT, &
            mesh_deformedCoordsLinear, &
            mesh_volumeMismatch, &
            mesh_shapeMismatch
#endif

 private :: &
#ifdef Spectral
            mesh_spectral_getGrid, &
            mesh_spectral_getSize, &
            mesh_spectral_getHomogenization, &
            mesh_spectral_count_nodesAndElements, &
            mesh_spectral_count_cpElements, &
            mesh_spectral_map_elements, &
            mesh_spectral_map_nodes, &
            mesh_spectral_count_cpSizes, &
            mesh_spectral_build_nodes, &
            mesh_spectral_build_elements, &
#endif 
#ifdef Marc
            mesh_marc_get_tableStyles, &
            mesh_marc_count_nodesAndElements, &
            mesh_marc_count_elementSets, &
            mesh_marc_map_elementSets, &
            mesh_marc_count_cpElements, &
            mesh_marc_map_Elements, &
            mesh_marc_map_nodes, &    
            mesh_marc_build_nodes, &
            mesh_marc_count_cpSizes, &
            mesh_marc_build_elements, &
#endif 
#ifdef Abaqus
            mesh_abaqus_count_nodesAndElements, &
            mesh_abaqus_count_elementSets, &
            mesh_abaqus_count_materials, &
            mesh_abaqus_map_elementSets, &
            mesh_abaqus_map_materials, &
            mesh_abaqus_count_cpElements, &
            mesh_abaqus_map_elements, &
            mesh_abaqus_map_nodes, &
            mesh_abaqus_build_nodes, &
            mesh_abaqus_count_cpSizes, &
            mesh_abaqus_build_elements, &
#endif 
            mesh_get_damaskOptions, &
            mesh_build_ipAreas, &
            mesh_build_nodeTwins, &
            mesh_build_sharedElems, &
            mesh_build_ipNeighborhood, &
            mesh_tell_statistics, &
            FE_mapElemtype, &
            mesh_faceMatch, &
            mesh_build_FEdata
contains


!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

 use DAMASK_interface
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_timeStamp, &
#ifdef Abaqus
   IO_abaqus_hasNoPart, &
#endif
#ifdef Spectral
   IO_open_file
 use numerics, only: &
   divergence_correction
#else
   IO_open_InputFile
#endif

 use FEsolving, only: &
   parallelExecution, &
   FEsolving_execElem, &
   FEsolving_execIP, &
   calcMode, &
   lastMode, &
   modelName
 
 implicit none
 integer(pInt), parameter :: fileUnit = 222_pInt
 integer(pInt), intent(in) :: el, ip
 integer(pInt) :: j
 
 write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 if (allocated(mesh_mapFEtoCPelem))  deallocate(mesh_mapFEtoCPelem)
 if (allocated(mesh_mapFEtoCPnode))  deallocate(mesh_mapFEtoCPnode)
 if (allocated(mesh_node0))          deallocate(mesh_node0)
 if (allocated(mesh_node))           deallocate(mesh_node)
 if (allocated(mesh_element))        deallocate(mesh_element)
 if (allocated(mesh_subNodeCoord))   deallocate(mesh_subNodeCoord)
 if (allocated(mesh_ipCoordinates))  deallocate(mesh_ipCoordinates)
 if (allocated(mesh_ipArea))         deallocate(mesh_ipArea)
 if (allocated(mesh_ipAreaNormal))   deallocate(mesh_ipAreaNormal)
 if (allocated(mesh_sharedElem))     deallocate(mesh_sharedElem)
 if (allocated(mesh_ipNeighborhood)) deallocate(mesh_ipNeighborhood)
 if (allocated(mesh_ipVolume))       deallocate(mesh_ipVolume)
 if (allocated(mesh_nodeTwins))      deallocate(mesh_nodeTwins)
 if (allocated(FE_nodesAtIP))        deallocate(FE_nodesAtIP)
 if (allocated(FE_ipNeighbor))       deallocate(FE_ipNeighbor)
 if (allocated(FE_subNodeParent))    deallocate(FE_subNodeParent)
 if (allocated(FE_subNodeOnIPFace))  deallocate(FE_subNodeOnIPFace)
 call mesh_build_FEdata                                                                             ! get properties of the different types of elements
#ifdef Spectral
 call IO_open_file(fileUnit,geometryFile)                                                           ! parse info from geometry file...
 res = mesh_spectral_getGrid(fileUnit)
 res1_red = res(1)/2_pInt + 1_pInt
 wgt = 1.0/real(product(res),pReal)
 geomdim = mesh_spectral_getSize(fileUnit)
 homog = mesh_spectral_getHomogenization(fileUnit)
 
!--------------------------------------------------------------------------------------------------
! scale dimension to calculate either uncorrected, dimension-independent, or dimension- and reso-
! lution-independent divergence
 if (divergence_correction == 1_pInt) then
   do j = 1_pInt, 3_pInt
    if (j /= minloc(geomdim,1) .and. j /= maxloc(geomdim,1)) scaledDim = geomdim/geomdim(j)
   enddo
 elseif (divergence_correction == 2_pInt) then
   do j = 1_pInt, 3_pInt
    if (j /= minloc(geomdim/res,1) .and. j /= maxloc(geomdim/res,1)) scaledDim = geomdim/geomdim(j)*res(j)
   enddo
 else
   scaledDim = geomdim
 endif
 write(6,'(a,3(i12  ))') ' grid     a b c: ', res
 write(6,'(a,3(f12.5))') ' size     x y z: ', geomdim
 write(6,'(a,i5,/)')     ' homogenization: ', homog
 call mesh_spectral_count_nodesAndElements
 call mesh_spectral_count_cpElements
 call mesh_spectral_map_elements
 call mesh_spectral_map_nodes
 call mesh_spectral_count_cpSizes
 call mesh_spectral_build_nodes
 call mesh_spectral_build_elements(fileUnit)
#endif
#ifdef Marc
 call IO_open_inputFile(fileUnit,modelName)                                                         ! parse info from input file...
 call mesh_marc_get_tableStyles(fileUnit)
 call mesh_marc_count_nodesAndElements(fileUnit)
 call mesh_marc_count_elementSets(fileUnit)
 call mesh_marc_map_elementSets(fileUnit)
 call mesh_marc_count_cpElements(fileUnit)
 call mesh_marc_map_elements(fileUnit)
 call mesh_marc_map_nodes(fileUnit)
 call mesh_marc_build_nodes(fileUnit)
 call mesh_marc_count_cpSizes(fileunit)
 call mesh_marc_build_elements(fileUnit)
#endif
#ifdef Abaqus
 call IO_open_inputFile(fileUnit,modelName)                                                         ! parse info from input file...
 noPart = IO_abaqus_hasNoPart(fileUnit)
 call mesh_abaqus_count_nodesAndElements(fileUnit)
 call mesh_abaqus_count_elementSets(fileUnit)
 call mesh_abaqus_count_materials(fileUnit)
 call mesh_abaqus_map_elementSets(fileUnit)
 call mesh_abaqus_map_materials(fileUnit)
 call mesh_abaqus_count_cpElements(fileUnit)
 call mesh_abaqus_map_elements(fileUnit)
 call mesh_abaqus_map_nodes(fileUnit)
 call mesh_abaqus_build_nodes(fileUnit)
 call mesh_abaqus_count_cpSizes(fileunit)
 call mesh_abaqus_build_elements(fileUnit)
#endif

 call mesh_get_damaskOptions(fileUnit)
 close (fileUnit)
 
 call mesh_build_subNodeCoords
 call mesh_build_ipCoordinates
 call mesh_build_ipVolumes
 call mesh_build_ipAreas
 call mesh_build_nodeTwins
 call mesh_build_sharedElems
 call mesh_build_ipNeighborhood
 call mesh_tell_statistics

 parallelExecution = (parallelExecution .and. (mesh_Nelems == mesh_NcpElems))                       ! plus potential killer from non-local constitutive
 
 FEsolving_execElem = [ 1_pInt,mesh_NcpElems ]
 if (allocated(FEsolving_execIP)) deallocate(FEsolving_execIP)
 allocate(FEsolving_execIP(2_pInt,mesh_NcpElems)); FEsolving_execIP = 1_pInt
 forall (j = 1_pInt:mesh_NcpElems) FEsolving_execIP(2,j) = FE_Nips(FE_geomtype(mesh_element(2,j)))
 
 if (allocated(calcMode)) deallocate(calcMode)
 allocate(calcMode(mesh_maxNips,mesh_NcpElems))
 calcMode = .false.                                                                                 ! pretend to have collected what first call is asking (F = I)
 calcMode(ip,mesh_FEasCP('elem',el)) = .true.                                                       ! first ip,el needs to be already pingponged to "calc"
 lastMode = .true.                                                                                  ! and its mode is already known...

end subroutine mesh_init

!--------------------------------------------------------------------------------------------------
!> @brief Gives the FE to CP ID mapping by binary search through lookup array
!! valid questions (what) are 'elem', 'node'
!--------------------------------------------------------------------------------------------------
integer(pInt) function mesh_FEasCP(what,myID)
 use IO, only: &
   IO_lc

 implicit none
 character(len=*), intent(in) :: what
 integer(pInt),    intent(in) :: myID
 
 integer(pInt), dimension(:,:), pointer :: lookupMap
 integer(pInt) :: lower,upper,center
 
 mesh_FEasCP = 0_pInt
 select case(IO_lc(what(1:4)))
   case('elem')
     lookupMap => mesh_mapFEtoCPelem
   case('node')
     lookupMap => mesh_mapFEtoCPnode
   case default
     return
 endselect
 
 lower = 1_pInt
 upper = int(size(lookupMap,2_pInt),pInt)
 
 if (lookupMap(1_pInt,lower) == myID) then                                                          ! check at bounds QUESTION is it valid to extend bounds by 1 and just do binary search w/o init check at bounds?
   mesh_FEasCP = lookupMap(2_pInt,lower)
   return
 elseif (lookupMap(1_pInt,upper) == myID) then
   mesh_FEasCP = lookupMap(2_pInt,upper)
   return
 endif
 
 do while (upper-lower > 1_pInt)                                                                    ! binary search in between bounds
   center = (lower+upper)/2_pInt
   if (lookupMap(1_pInt,center) < myID) then
     lower = center
   elseif (lookupMap(1_pInt,center) > myID) then
     upper = center
   else
     mesh_FEasCP = lookupMap(2_pInt,center)
     exit
   endif
 enddo
 
end function mesh_FEasCP


!--------------------------------------------------------------------------------------------------
!> @brief Assigns coordinates for subnodes in each CP element.
!! Allocates global array 'mesh_subNodeCoord'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_subNodeCoords
 
 implicit none
 integer(pInt) e,t,n,p,Nparents
 real(pReal), dimension(3,mesh_maxNnodes+mesh_maxNsubNodes) :: mySubNodeCoord
 
 if (.not. allocated(mesh_subNodeCoord)) then
   allocate(mesh_subNodeCoord(3,mesh_maxNnodes+mesh_maxNsubNodes,mesh_NcpElems))
 endif
 mesh_subNodeCoord = 0.0_pReal
 
 !$OMP PARALLEL DO PRIVATE(mySubNodeCoord,t,Nparents)
 do e = 1_pInt,mesh_NcpElems                                                                        ! loop over cpElems
   mySubNodeCoord = 0.0_pReal
   t = FE_geomtype(mesh_element(2,e))                                                               ! get elemGeomType
   do n = 1_pInt,FE_Nnodes(t)
     mySubNodeCoord(1:3,n) = mesh_node(1:3,mesh_FEasCP('node',mesh_element(4_pInt+n,e)))            ! loop over nodes of this element type
   enddo
   do n = 1_pInt,FE_NsubNodes(t)                                                                    ! now for the true subnodes
     Nparents = count(FE_subNodeParent(1_pInt:FE_Nips(t),n,t) > 0_pInt)
     do p = 1_pInt,Nparents                                                                         ! loop through present parent nodes
       mySubNodeCoord(1:3,FE_Nnodes(t)+n) &
                = mySubNodeCoord(1:3,FE_Nnodes(t)+n) &
                + mesh_node(1:3,mesh_FEasCP('node',mesh_element(4_pInt+FE_subNodeParent(p,n,t),e))) ! add up parents
     enddo
     mySubNodeCoord(1:3,n+FE_Nnodes(t)) = mySubNodeCoord(1:3,n+FE_Nnodes(t)) / real(Nparents,pReal)
   enddo
   mesh_subNodeCoord(1:3,1:mesh_maxNnodes+mesh_maxNsubNodes,e) = mySubNodeCoord
 enddo 
 !$OMP END PARALLEL DO
 
end subroutine mesh_build_subNodeCoords


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume. Allocates global array 'mesh_ipVolume'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipVolumes
 use math, only: &
   math_volTetrahedron, &
   math_areaTriangle
 
 implicit none
 integer(pInt) ::                                e,f,t,i,j,n, &
                                                 NmySubNodes, &              !< number of subnodes already found for this (2d) ip cell
                                                 Ntriangles                  !< each interface is made up of this many triangles
 integer(pInt), dimension(FE_maxNipFaceNodes) :: mySubNodes                  !< list of subnodes that constitute this (2d) ip cell
 real(pReal)                                  :: myVolume
 real(pReal), dimension(3,FE_maxNipFaceNodes) :: nPos                        !< coordinates of nodes on IP face
 real(pReal), dimension(FE_maxNipFaceNodes-2_pInt,FE_maxNipFaceNodes) :: volume

 if (.not. allocated(mesh_ipVolume)) then
   allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems))
 endif
 
 mesh_ipVolume = 0.0_pReal 
!$OMP PARALLEL DO PRIVATE(t,Ntriangles,nPos,volume,mySubNodes,NmySubNodes,myVolume)
 do e = 1_pInt,mesh_NcpElems                                                                        ! loop over cpElems
   t = FE_geomtype(mesh_element(2,e))                                                               ! get elemGeomType
   select case (FE_dimension(t))                                                                    ! treat (1D), 2D, and 3D element types differently

     case (2_pInt)
       Ntriangles = FE_NipNeighbors(t) - 2_pInt
       do i = 1_pInt,FE_Nips(t)                                                                     ! loop over IPs of elem
         nPos = 0.0_pReal
         volume = 0.0_pReal
         mySubNodes = 0_pInt
         NmySubNodes = 0_pInt
faceLoop:do f = 1_pInt,FE_NipNeighbors(t)                                                           ! loop over interfaces of IP and add tetrahedra which connect to CoG
           do n = 1_pInt,FE_NipFaceNodes(t)
             if (all(mySubNodes /= FE_subNodeOnIPFace(n,f,i,t))) then
               NmySubNodes = NmySubNodes + 1_pInt
               mySubNodes(NmySubNodes) = FE_subNodeOnIPFace(n,f,i,t)
               nPos(1:3,NmySubNodes) = mesh_subNodeCoord(1:3,mySubNodes(NmySubNodes),e)
             elseif (NmySubNodes == Fe_maxNipFaceNodes) then
               exit faceLoop
             endif
           enddo
         enddo faceLoop
         forall (n = 1_pInt:FE_NipNeighbors(t), j = 1_pInt:Ntriangles) &                            ! start at each interface node and build valid triangles to cover interface
           volume(j,n) = math_areaTriangle(nPos(1:3,n), &
                                           nPos(1:3,1_pInt+mod(n+j-1_pInt,FE_NipNeighbors(t))), &
                                           nPos(1:3,1_pInt+mod(n+j       ,FE_NipNeighbors(t))))     ! assuming element depth of 1
         mesh_ipVolume(i,e) = sum(volume) / FE_NipFaceNodes(t)                                      ! renormalize with interfaceNodeNum due to loop over them
       enddo

     case (3_pInt)
       Ntriangles = FE_NipFaceNodes(t) - 2_pInt
       do i = 1_pInt,FE_Nips(t)                                                                     ! loop over IPs of elem
         myVolume = 0.0_pReal
         do f = 1_pInt,FE_NipNeighbors(t)                                                           ! loop over interfaces of IP and add tetrahedra which connect to CoG
           nPos = 0.0_pReal
           volume = 0.0_pReal
           forall (n = 1_pInt:FE_NipFaceNodes(t)) &
             nPos(:,n) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
           forall (n = 1_pInt:FE_NipFaceNodes(t), j = 1_pInt:Ntriangles) &                          ! start at each interface node and build valid triangles to cover interface
             volume(j,n) = math_volTetrahedron(nPos(:,n), &                                         ! calc volume of respective tetrahedron to CoG
                                               nPos(:,1_pInt+mod(n-1_pInt +j       ,FE_NipFaceNodes(t))),& ! start at offset j
                                               nPos(:,1_pInt+mod(n-1_pInt +j+1_pInt,FE_NipFaceNodes(t))),& ! and take j's neighbor
                                               mesh_cellCenterCoordinates(i,e))
           myVolume = myVolume + sum(volume)                                                        ! add contribution from this interface
         enddo
         mesh_ipVolume(i,e) = myVolume / FE_NipFaceNodes(t)                                         ! renormalize with interfaceNodeNum due to loop over them
       enddo
   endselect
 enddo
!$OMP END PARALLEL DO

end subroutine mesh_build_ipVolumes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP Coordinates. Allocates global array 'mesh_ipCoordinates'
! Called by all solvers in mesh_init in order to initialize the ip coordinates.
! Later on the current ip coordinates are directly prvided by the spectral solver and by Abaqus,
! so no need to use this subroutine anymore; Marc however only provides nodal displacements,
! so in this case the ip coordinates are always calculated on the basis of this subroutine.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR THE MOMENT THIS SUBROUTINE ACTUALLY CALCULATES THE CELL CELLENTER AND NOT THE IP COORDINATES,
! AS THE IP IS NOT (ALWAYS) LOCATED IN THE CENTER OF THE IP VOLUME. 
! HAS TO BE CHANGED IN A LATER VERSION. 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipCoordinates
 use prec, only: &
   tol_gravityNodePos
 
 implicit none
 integer(pInt) :: e,f,t,i,j,k,n
 logical, dimension(mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNode                                ! flagList to find subnodes determining center of grav
 real(pReal), dimension(3,mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNodePos                       ! coordinates of subnodes determining center of grav

 if (.not. allocated(mesh_ipCoordinates)) allocate(mesh_ipCoordinates(3,mesh_maxNips,mesh_NcpElems))
 
 !$OMP PARALLEL DO PRIVATE(t,gravityNode,gravityNodePos)
 do e = 1_pInt,mesh_NcpElems                                                                       ! loop over cpElems
   t = FE_geomtype(mesh_element(2,e))                                                              ! get elemGeomType
   do i = 1_pInt,FE_Nips(t)                                                                        ! loop over IPs of elem
     gravityNode = .false.                                                                         ! reset flagList
     gravityNodePos = 0.0_pReal                                                                    ! reset coordinates
     do f = 1_pInt,FE_NipNeighbors(t)                                                              ! loop over interfaces of IP
       do n = 1_pInt,FE_NipFaceNodes(t)                                                            ! loop over nodes on interface
         gravityNode(FE_subNodeOnIPFace(n,f,i,t)) = .true.
         gravityNodePos(:,FE_subNodeOnIPFace(n,f,i,t)) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       enddo
     enddo
     
     do j = 1_pInt,mesh_maxNnodes+mesh_maxNsubNodes-1_pInt                                         ! walk through entire flagList except last
       if (gravityNode(j)) then                                                                    ! valid node index
         do k = j+1_pInt,mesh_maxNnodes+mesh_maxNsubNodes                                          ! walk through remainder of list
           if (gravityNode(k) .and. all(abs(gravityNodePos(:,j) - gravityNodePos(:,k)) < tol_gravityNodePos)) then   ! found duplicate
             gravityNode(j) = .false.                                                              ! delete first instance
             gravityNodePos(:,j) = 0.0_pReal
             exit                                                                                  ! continue with next suspect
           endif
         enddo
       endif
     enddo
     mesh_ipCoordinates(:,i,e) = sum(gravityNodePos,2)/real(count(gravityNode),pReal)
   enddo
 enddo
 !$OMP END PARALLEL DO

end subroutine mesh_build_ipCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief Calculates cell center coordinates.
!--------------------------------------------------------------------------------------------------
pure function mesh_cellCenterCoordinates(i,e)
use prec, only: &
  tol_gravityNodePos
 
implicit none
!*** input variables
integer(pInt), intent(in) :: e, &                                                                  ! element number
                             i                                                                     ! integration point number

!*** output variables
real(pReal), dimension(3) :: mesh_cellCenterCoordinates                                            ! x,y,z coordinates of the cell center of the requested IP cell

!*** local variables
integer(pInt) :: f,t,j,k,n
logical, dimension(mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNode                                ! flagList to find subnodes determining center of grav
real(pReal), dimension(3,mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNodePos                       ! coordinates of subnodes determining center of grav
 

t = FE_geomtype(mesh_element(2,e))                                                                 ! get elemGeomType
gravityNode = .false.                                                                              ! reset flagList
gravityNodePos = 0.0_pReal                                                                         ! reset coordinates
do f = 1_pInt,FE_NipNeighbors(t)                                                                   ! loop over interfaces of IP
  do n = 1_pInt,FE_NipFaceNodes(t)                                                                    ! loop over nodes on interface
    gravityNode(FE_subNodeOnIPFace(n,f,i,t)) = .true.
    gravityNodePos(:,FE_subNodeOnIPFace(n,f,i,t)) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
  enddo
enddo
do j = 1_pInt,mesh_maxNnodes+mesh_maxNsubNodes-1_pInt                                              ! walk through entire flagList except last
  if (gravityNode(j)) then                                                                         ! valid node index
    do k = j+1_pInt,mesh_maxNnodes+mesh_maxNsubNodes                                               ! walk through remainder of list
      if (gravityNode(k) .and. all(abs(gravityNodePos(:,j) - gravityNodePos(:,k)) < tol_gravityNodePos)) then   ! found duplicate
        gravityNode(j) = .false.                                                                   ! delete first instance
        gravityNodePos(:,j) = 0.0_pReal
        exit                                                                                       ! continue with next suspect
      endif
    enddo
  endif
enddo
mesh_cellCenterCoordinates = sum(gravityNodePos,2)/real(count(gravityNode),pReal)

endfunction mesh_cellCenterCoordinates


#ifdef Spectral
!--------------------------------------------------------------------------------------------------
!> @brief Reads grid information from geometry file. If fileUnit is given, 
!! assumes an opened file, otherwise tries to open the one specified in geometryFile
!--------------------------------------------------------------------------------------------------
function mesh_spectral_getGrid(fileUnit)
 use IO, only: &
   IO_checkAndRewind, &
   IO_open_file, &
   IO_stringPos, &
   IO_lc, &
   IO_stringValue, &
   IO_intValue, &
   IO_floatValue, &
   IO_error
 use DAMASK_interface, only: &
   geometryFile
  
 implicit none
 integer(pInt), dimension(3)                      :: mesh_spectral_getGrid
 integer(pInt), intent(in), optional              :: fileUnit
 integer(pInt), dimension(1_pInt + 7_pInt*2_pInt) :: positions                                     ! for a,b,c + 3 values + keyword

 integer(pInt)                                    :: headerLength = 0_pInt
 character(len=1024) :: line, &
                        keyword
 integer(pInt) :: i, j, myUnit
 logical :: gotGrid = .false.
 
 mesh_spectral_getGrid = -1_pInt
 if(.not. present(fileUnit)) then
   myUnit = 289_pInt
   call IO_open_file(myUnit,trim(geometryFile))
 else
   myUnit = fileUnit
 endif
 
 call IO_checkAndRewind(myUnit)

 read(myUnit,'(a1024)') line
 positions = IO_stringPos(line,7_pInt)
 keyword = IO_lc(IO_StringValue(line,positions,2_pInt,.true.))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,positions,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=841_pInt, ext_msg='mesh_spectral_getGrid')
 endif
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   positions = IO_stringPos(line,7_pInt)             
   select case ( IO_lc(IO_StringValue(line,positions,1_pInt,.true.)) )
     case ('resolution','grid')
       gotGrid = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,positions,j)))
           case('a')
              mesh_spectral_getGrid(1) = IO_intValue(line,positions,j+1_pInt)
           case('b')
              mesh_spectral_getGrid(2) = IO_intValue(line,positions,j+1_pInt)
           case('c')
              mesh_spectral_getGrid(3) = IO_intValue(line,positions,j+1_pInt)
         end select
       enddo
   end select
 enddo
 
 if(.not. present(fileUnit)) close(myUnit)
 
 if (.not. gotGrid) &
   call IO_error(error_ID = 845_pInt, ext_msg='grid')
 if(any(mesh_spectral_getGrid < 1_pInt)) &
   call IO_error(error_ID = 843_pInt, ext_msg='mesh_spectral_getGrid')

end function mesh_spectral_getGrid


!--------------------------------------------------------------------------------------------------
!> @brief Reads size information from geometry file. If fileUnit is given, 
!! assumes an opened file, otherwise tries to open the one specified in geometryFile
!--------------------------------------------------------------------------------------------------
function mesh_spectral_getSize(fileUnit)
 use IO, only: &
   IO_checkAndRewind, &
   IO_open_file, &
   IO_stringPos, &
   IO_lc, &
   IO_stringValue, &
   IO_intValue, &
   IO_floatValue, &
   IO_error
 use DAMASK_interface, only: &
   geometryFile
  
 implicit none
 real(pReal), dimension(3)                        :: mesh_spectral_getSize
 integer(pInt), intent(in), optional              :: fileUnit
 integer(pInt), dimension(1_pInt + 7_pInt*2_pInt) :: positions                                      ! for x,y,z + 3 values + keyword
 integer(pInt)                                    :: headerLength = 0_pInt
 character(len=1024) :: line, &
                        keyword
 integer(pInt) :: i, j, myUnit 
 logical :: gotSize = .false.
 
 mesh_spectral_getSize = -1.0_pReal
 if(.not. present(fileUnit)) then
   myUnit = 289_pInt
   call IO_open_file(myUnit,trim(geometryFile))
 else
   myUnit = fileUnit
 endif
 
 call IO_checkAndRewind(myUnit)

 read(myUnit,'(a1024)') line
 positions = IO_stringPos(line,7_pInt)
 keyword = IO_lc(IO_StringValue(line,positions,2_pInt,.true.))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,positions,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=841_pInt, ext_msg='mesh_spectral_getSize')
 endif
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   positions = IO_stringPos(line,7_pInt)             
   select case ( IO_lc(IO_StringValue(line,positions,1,.true.)) )
     case ('dimension', 'size')
       gotSize = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,positions,j)))
           case('x')
              mesh_spectral_getSize(1) = IO_floatValue(line,positions,j+1_pInt)
           case('y')
              mesh_spectral_getSize(2) = IO_floatValue(line,positions,j+1_pInt)
           case('z')
              mesh_spectral_getSize(3) = IO_floatValue(line,positions,j+1_pInt)
         end select
       enddo
   end select
 enddo
 
 if(.not. present(fileUnit)) close(myUnit)

 if (.not. gotSize) &
   call IO_error(error_ID = 845_pInt, ext_msg='size')
 if (any(mesh_spectral_getSize<=0.0_pReal)) &
   call IO_error(error_ID = 844_pInt, ext_msg='mesh_spectral_getSize')

end function mesh_spectral_getSize


!--------------------------------------------------------------------------------------------------
!> @brief Reads homogenization information from geometry file. If fileUnit is given, 
!! assumes an opened file, otherwise tries to open the one specified in geometryFile
!--------------------------------------------------------------------------------------------------
integer(pInt) function mesh_spectral_getHomogenization(fileUnit)
 use IO, only: &
   IO_checkAndRewind, &
   IO_open_file, &
   IO_stringPos, &
   IO_lc, &
   IO_stringValue, &
   IO_intValue, &
   IO_error
 use DAMASK_interface, only: &
   geometryFile
  
 implicit none
 integer(pInt), intent(in), optional              :: fileUnit
 integer(pInt), dimension(1_pInt + 7_pInt*2_pInt) :: positions                                      ! for a, b,  c + 3 values + keyword
 integer(pInt)                                    :: headerLength = 0_pInt
 character(len=1024) :: line, &
                        keyword
 integer(pInt) :: i, myUnit
 logical :: gotHomogenization = .false.
 
 mesh_spectral_getHomogenization = -1_pInt
 if(.not. present(fileUnit)) then
   myUnit = 289_pInt
   call IO_open_file(myUnit,trim(geometryFile))
 else
   myUnit = fileUnit
 endif
 
 call IO_checkAndRewind(myUnit)

 read(myUnit,'(a1024)') line
 positions = IO_stringPos(line,7_pInt)
 keyword = IO_lc(IO_StringValue(line,positions,2_pInt,.true.))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,positions,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=841_pInt, ext_msg='mesh_spectral_getHomogenization')
 endif
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   positions = IO_stringPos(line,7_pInt)             
   select case ( IO_lc(IO_StringValue(line,positions,1,.true.)) )
     case ('homogenization')
       gotHomogenization = .true.
       mesh_spectral_getHomogenization = IO_intValue(line,positions,2_pInt)
   end select
 enddo
 
 if(.not. present(fileUnit)) close(myUnit)
 
 if (.not. gotHomogenization ) &
   call IO_error(error_ID = 845_pInt, ext_msg='homogenization')
 if (mesh_spectral_getHomogenization<1_pInt) &
   call IO_error(error_ID = 842_pInt, ext_msg='mesh_spectral_getHomogenization')
   
end function mesh_spectral_getHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of nodes and elements in mesh and stores them in
!! 'mesh_Nelems' and 'mesh_Nnodes'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_count_nodesAndElements()
 
 implicit none
 mesh_Nelems = product(res)
 mesh_Nnodes = product(res+1_pInt)

end subroutine mesh_spectral_count_nodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of CP elements in mesh and stores them in 'mesh_NcpElems'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_count_cpElements

 implicit none

 mesh_NcpElems = mesh_Nelems
 
end subroutine mesh_spectral_count_cpElements


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPelem'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_map_elements

 implicit none
 integer(pInt) :: i

 allocate (mesh_mapFEtoCPelem(2_pInt,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

 forall (i = 1_pInt:mesh_NcpElems) &
   mesh_mapFEtoCPelem(1:2,i) = i

end subroutine mesh_spectral_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPnode'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_map_nodes

 implicit none
 integer(pInt) :: i

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

 forall (i = 1_pInt:mesh_Nnodes) &
   mesh_mapFEtoCPnode(1:2,i) = i
 
end subroutine mesh_spectral_map_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets maximum count of nodes, IPs, IP neighbors, and subNodes among cpElements.
!! Allocates global arrays 'mesh_maxNnodes', 'mesh_maxNips', mesh_maxNipNeighbors', 
!! and mesh_maxNsubNodes
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_count_cpSizes
 
 implicit none
 integer(pInt) :: t
 
 t = FE_geomtype(FE_mapElemtype('C3D8R'))                                                             ! fake 3D hexahedral 8 node 1 IP element

 mesh_maxNnodes =       FE_Nnodes(t)
 mesh_maxNips =         FE_Nips(t)
 mesh_maxNipNeighbors = FE_NipNeighbors(t)
 mesh_maxNsubNodes =    FE_NsubNodes(t)

end subroutine mesh_spectral_count_cpSizes


!--------------------------------------------------------------------------------------------------
!> @brief Store x,y,z coordinates of all nodes in mesh.
!! Allocates global arrays 'mesh_node0' and 'mesh_node'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_nodes()
 use numerics, &
   only: numerics_unitlength

 implicit none
 integer(pInt) :: n

 allocate ( mesh_node0 (3,mesh_Nnodes) ); mesh_node0 = 0.0_pReal
 allocate ( mesh_node  (3,mesh_Nnodes) ); mesh_node  = 0.0_pReal
 
 forall (n = 0_pInt:mesh_Nnodes-1_pInt)
   mesh_node0(1,n+1_pInt) = numerics_unitlength * &
           geomdim(1) * real(mod(n,(res(1)+1_pInt) ),pReal) &
                                                   / real(res(1),pReal)
   mesh_node0(2,n+1_pInt) = numerics_unitlength * &
           geomdim(2) * real(mod(n/(res(1)+1_pInt),(res(2)+1_pInt)),pReal) &
                                                   / real(res(2),pReal)
   mesh_node0(3,n+1_pInt) = numerics_unitlength * &
           geomdim(3) * real(mod(n/(res(1)+1_pInt)/(res(2)+1_pInt),(res(3)+1_pInt)),pReal) &
                                                   / real(res(3),pReal)
 end forall 

 mesh_node = mesh_node0                                                                             ! why?

end subroutine mesh_spectral_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, material, texture, and node list per element.
!! Allocates global array 'mesh_element'
!--------------------------------------------------------------------------------------------------
subroutine mesh_spectral_build_elements(myUnit)

 use IO, only: &
   IO_checkAndRewind, &
   IO_lc, &
   IO_stringValue, &
   IO_stringPos, &
   IO_error, &
   IO_continuousIntValues, &
   IO_intValue, &
   IO_countContinuousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), dimension (1_pInt+7_pInt*2_pInt) :: myPos
 integer(pInt) :: e, i, headerLength = 0_pInt, maxIntCount
 integer(pInt), dimension(:), allocatable :: microstructures
 integer(pInt), dimension(1,1) :: dummySet = 0_pInt
 character(len=65536) :: line,keyword
 character(len=64), dimension(1) :: dummyName = ''

 call IO_checkAndRewind(myUnit)

 read(myUnit,'(a65536)') line
 myPos = IO_stringPos(line,7_pInt)
 keyword = IO_lc(IO_StringValue(line,myPos,2_pInt,.true.))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,myPos,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=841_pInt, ext_msg='mesh_spectral_build_elements')
 endif
 
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a65536)') line
 enddo

 maxIntCount = 0_pInt
 i = 1_pInt

 do while (i > 0_pInt)
   i = IO_countContinuousIntValues(myUnit)
   maxIntCount = max(maxIntCount, i)
 enddo

 rewind (myUnit)
 do i=1_pInt,headerLength                                                                           ! skip header
   read(myUnit,'(a65536)') line
 enddo

 allocate (mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems)); mesh_element = 0_pInt
 allocate (microstructures (1_pInt+maxIntCount));               microstructures = 2_pInt
 
 e = 0_pInt
 do while (e < mesh_NcpElems .and. microstructures(1) > 0_pInt)                                     ! fill expected number of elements, stop at end of data (or blank line!)
   microstructures = IO_continuousIntValues(myUnit,maxIntCount,dummyName,dummySet,0_pInt)           ! get affected elements
   do i = 1_pInt,microstructures(1_pInt)
     e = e+1_pInt                                                                                   ! valid element entry
     mesh_element( 1,e) = e                                                                         ! FE id
     mesh_element( 2,e) = FE_mapElemtype('C3D8R')                                                   ! elem type
     mesh_element( 3,e) = homog                                                                     ! homogenization
     mesh_element( 4,e) = microstructures(1_pInt+i)                                                 ! microstructure
     mesh_element( 5,e) = e + (e-1_pInt)/res(1) + &
                                       ((e-1_pInt)/(res(1)*res(2)))*(res(1)+1_pInt)                 ! base node
     mesh_element( 6,e) = mesh_element(5,e) + 1_pInt
     mesh_element( 7,e) = mesh_element(5,e) + res(1) + 2_pInt
     mesh_element( 8,e) = mesh_element(5,e) + res(1) + 1_pInt
     mesh_element( 9,e) = mesh_element(5,e) +(res(1) + 1_pInt) * (res(2) + 1_pInt)                  ! second floor base node
     mesh_element(10,e) = mesh_element(9,e) + 1_pInt
     mesh_element(11,e) = mesh_element(9,e) + res(1) + 2_pInt
     mesh_element(12,e) = mesh_element(9,e) + res(1) + 1_pInt
     mesh_maxValStateVar(1) = max(mesh_maxValStateVar(1),mesh_element(3,e))                         !needed for statistics
     mesh_maxValStateVar(2) = max(mesh_maxValStateVar(2),mesh_element(4,e))              
   enddo
 enddo

 deallocate(microstructures)
 if (e /= mesh_NcpElems) call IO_error(880_pInt,e)                                                  !@ToDo does that make sense?

end subroutine mesh_spectral_build_elements


!--------------------------------------------------------------------------------------------------
!> @brief Performes a regridding from saved restart information
!--------------------------------------------------------------------------------------------------
function mesh_regrid(adaptive,resNewInput,minRes)
 use prec, only: &
   pInt, &
   pReal
 use DAMASK_interface, only: &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   GeometryFile
 use IO, only: &
   IO_read_jobBinaryFile ,&
   IO_read_jobBinaryIntFile ,&
   IO_write_jobBinaryFile, &
   IO_write_jobBinaryIntFile, &
   IO_write_jobFile, &
   IO_error
 use numerics, only: &
   mySpectralSolver
 use math, only: &
   math_periodicNearestNeighbor, &
   math_mul33x3
 character(len=1024):: formatString, N_Digits
 logical, intent(in)                                    :: adaptive                                  ! if true, choose adaptive grid based on resNewInput, otherwise keep it constant
 integer(pInt), dimension(3), optional, intent(in)      :: resNewInput                               ! f2py cannot handle optional arguments correctly (they are always present)
 integer(pInt), dimension(3), optional, intent(in)      :: minRes
 integer(pInt), dimension(3)                            :: mesh_regrid, ratio
 integer(pInt), dimension(3,2)                          :: possibleResNew
 integer(pInt):: maxsize, i, j, k, ielem, NpointsNew, spatialDim
 integer(pInt), dimension(3)                            :: resNew
 integer(pInt), dimension(:),            allocatable    :: indices
 real(pReal),   dimension(3)                            :: geomdimNew                                              
 real(pReal),   dimension(3,3)                          :: Favg, Favg_LastInc,       &
                                                           FavgNew, Favg_LastIncNew, &
                                                           deltaF, deltaF_lastInc
 real(pReal),   dimension(:,:),           allocatable :: & 
   coordinates,   coordinatesNew  
 real(pReal),   dimension(:,:,:),         allocatable :: & 
   stateHomog
 real(pReal),   dimension (:,:,:,:),      allocatable :: &
   spectralF9,   spectralF9New, &
   Tstar,           TstarNew, &
   stateConst 
 real(pReal),   dimension(:,:,:,:,:),     allocatable :: & 
   spectralF33,  spectralF33New, &
   F,                  FNew, &
   Fp,                FpNew, &
   Lp,                LpNew, &
   dcsdE,          dcsdENew, &
   F_lastIncNew
 real(pReal),   dimension (:,:,:,:,:,:,:), allocatable :: &
   dPdF,            dPdFNew

 integer(pInt), dimension(:,:), allocatable :: &
   sizeStateHomog
 integer(pInt), dimension(:,:,:), allocatable :: &
   material_phase, material_phaseNew, &
   sizeStateConst
 
 write(6,'(a)') 'Regridding geometry'
 if (adaptive) then
   write(6,'(a)') 'adaptive resolution determination'
   if (present(minRes)) then
     if (all(minRes /= -1_pInt)) &                                                                  !the f2py way to tell it is present
       write(6,'(a,3(i12))') ' given minimum resolution ', minRes
   endif
   if (present(resNewInput)) then
     if (any (resNewInput<1)) call IO_error(890_pInt, ext_msg = 'resNewInput')                      !the f2py way to tell it is not present
     write(6,'(a,3(i12))') ' target resolution ', resNewInput
   else
     call IO_error(890_pInt, ext_msg = 'resNewInput')
   endif
 endif
 
 allocate(coordinates(3,mesh_NcpElems))
 
!--------------------------------------------------------------------------------------------------
! read in deformation gradient to calculate coordinates, shape depend of selected solver
 select case(myspectralsolver)
   case('basic')
     allocate(spectralF33(3,3,res(1),res(2),res(3)))
     call IO_read_jobBinaryFile(777,'F',trim(getSolverJobName()),size(spectralF33))
     read (777,rec=1) spectralF33
     close (777)
     Favg = sum(sum(sum(spectralF33,dim=5),dim=4),dim=3) * wgt
     coordinates = reshape(mesh_deformedCoordsFFT(geomdim,spectralF33),[3,mesh_NcpElems])
   case('basicpetsc','al')
     allocate(spectralF9(9,res(1),res(2),res(3)))
     call IO_read_jobBinaryFile(777,'F',trim(getSolverJobName()),size(spectralF9))
     read (777,rec=1) spectralF9
     close (777)
     Favg = reshape(sum(sum(sum(spectralF9,dim=4),dim=3),dim=2) * wgt, [3,3])
     coordinates = reshape(mesh_deformedCoordsFFT(geomdim,reshape(spectralF9, &
                                      [3,3,res(1),res(2),res(3)])),[3,mesh_NcpElems])
  end select
  
!--------------------------------------------------------------------------------------------------
!  sanity check 2D/3D case
 if (res(3)== 1_pInt) then
   spatialDim = 2_pInt
   if (present (minRes)) then
     if (minRes(1) > 0_pInt .or. minRes(2) > 0_pInt) then
        if (minRes(3) /= 1_pInt .or. &
           mod(minRes(1),2_pInt) /= 0_pInt .or. &
           mod(minRes(2),2_pInt) /= 0_pInt)  call IO_error(890_pInt, ext_msg = '2D minRes')                       ! as f2py has problems with present, use pyf file for initialization to -1
   endif; endif
 else
   spatialDim = 3_pInt
   if (present (minRes)) then
     if (any(minRes > 0_pInt)) then
        if (mod(minRes(1),2_pInt) /= 0_pInt.or. &
            mod(minRes(2),2_pInt) /= 0_pInt .or. &
            mod(minRes(3),2_pInt) /= 0_pInt)  call IO_error(890_pInt, ext_msg = '3D minRes')                      ! as f2py has problems with present, use pyf file for initialization to -1
   endif; endif
 endif
 
!--------------------------------------------------------------------------------------------------
!  Automatic detection based on current geom
 geomdimNew =  math_mul33x3(Favg,geomdim)
 if (adaptive) then
   ratio = floor(real(resNewInput,pReal) * (geomdimNew/geomdim), pInt)
   
   possibleResNew = 1_pInt
   do i = 1_pInt, spatialDim
     if (mod(ratio(i),2) == 0_pInt) then
       possibleResNew(i,1:2) = [ratio(i),ratio(i) + 2_pInt]
     else
       possibleResNew(i,1:2) = [ratio(i)-1_pInt, ratio(i) + 1_pInt]
     endif
     if (.not.present(minRes)) then                              ! calling from fortran, optional argument not given
       possibleResNew = possibleResNew
     else                                                        ! optional argument is there
       if (any(minRes<1_pInt)) then
         possibleResNew = possibleResNew                     ! f2py calling, but without specification (or choosing invalid values), standard from pyf = -1
       else                                                        ! given useful values
         do k = 1_pInt,3_pInt; do j = 1_pInt,3_pInt
            possibleResNew(j,k) = max(possibleResNew(j,k), minRes(j))
         enddo; enddo
       endif
     endif
   enddo
   
   k = huge(1_pInt)
   do i = 0_pInt, 2_pInt**spatialDim - 1
      j = abs( possibleResNew(1,iand(i,1_pInt)/1_pInt + 1_pInt) &
             * possibleResNew(2,iand(i,2_pInt)/2_pInt + 1_pInt) &
             * possibleResNew(3,iand(i,4_pInt)/4_pInt + 1_pInt) &
             - resNewInput(1)*resNewInput(2)*resNewInput(3))
       
     if (j < k) then 
       k = j
       resNew =[ possibleResNew(1,iand(i,1_pInt)/1_pInt + 1_pInt), &
                 possibleResNew(2,iand(i,2_pInt)/2_pInt + 1_pInt), &
                 possibleResNew(3,iand(i,4_pInt)/4_pInt + 1_pInt) ] 
     endif
   enddo 
 else 
  resNew = res
 endif

 mesh_regrid = resNew
 NpointsNew = resNew(1)*resNew(2)*resNew(3)

!--------------------------------------------------------------------------------------------------
!  Calculate regular new coordinates
 allocate(coordinatesNew(3,NpointsNew))
 ielem = 0_pInt
 do k=1_pInt,resNew(3); do j=1_pInt, resNew(2); do i=1_pInt, resNew(1)
   ielem = ielem + 1_pInt
   coordinatesNew(1:3,ielem) = math_mul33x3(Favg,  geomdim/real(resNew,pReal)*real([i,j,k],pReal) &
                                                 - geomdim/real(2_pInt*resNew,pReal))
 enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
!  Nearest neighbour search
 allocate(indices(NpointsNew))
 indices =  math_periodicNearestNeighbor(geomdim, Favg, coordinatesNew, coordinates)
 deallocate(coordinates)

!--------------------------------------------------------------------------------------------------
!  write out indices periodic
 write(N_Digits, '(I16.16)') 1_pInt + int(log10(real(maxval(indices),pReal)))
 N_Digits = adjustl(N_Digits)
 formatString = '(I'//trim(N_Digits)//'.'//trim(N_Digits)//',a)'

 call IO_write_jobFile(777,'IDX')                                                                   ! make it a general open-write file
 write(777, '(A)') '1 header'
 write(777, '(A)') 'Numbered indices as per the large set'
 do i = 1_pInt, NpointsNew
   write(777,trim(formatString),advance='no') indices(i), ' '
   if(mod(i,resNew(1)) == 0_pInt) write(777,'(A)') ''
 enddo
 close(777)
 
 
!--------------------------------------------------------------------------------------------------
! calculate and write out indices non periodic
 do i = 1_pInt, NpointsNew
   indices(i) = indices(i) / 3_pInt**spatialDim +1_pInt                                             ! +1 b'coz index count starts from '0'
 enddo 
 write(N_Digits, '(I16.16)') 1_pInt + int(log10(real(maxval(indices),pReal)))
 N_Digits = adjustl(N_Digits)
 formatString = '(I'//trim(N_Digits)//'.'//trim(N_Digits)//',a)'

 call IO_write_jobFile(777,'idx')                                                                   ! make it a general open-write file
 write(777, '(A)') '1 header'
 write(777, '(A)') 'Numbered indices as per the small set'
 do i = 1_pInt, NpointsNew
   write(777,trim(formatString),advance='no') indices(i), ' '
   if(mod(i,resNew(1)) == 0_pInt) write(777,'(A)') ''
 enddo
 close(777)

!--------------------------------------------------------------------------------------------------
! write out new geom file
 write(N_Digits, '(I16.16)') 1_pInt+int(log10(real(maxval(mesh_element(4,1:mesh_NcpElems)),pReal)),pInt)
 N_Digits = adjustl(N_Digits)
 formatString = '(I'//trim(N_Digits)//'.'//trim(N_Digits)//',a)'
 open(777,file=trim(getSolverWorkingDirectoryName())//trim(GeometryFile),status='REPLACE')
 write(777, '(A)') '3 header'
 write(777, '(3(A, I8))') 'resolution  a ', resNew(1), '  b ', resNew(2), '  c ', resNew(3)
 write(777, '(3(A, g17.10))') 'dimension   x ', geomdim(1), '  y ', geomdim(2), '  z ', geomdim(3)
 write(777, '(A)') 'homogenization  1'
 do i = 1_pInt, NpointsNew
   write(777,trim(formatString),advance='no') mesh_element(4,indices(i)), ' '
   if(mod(i,resNew(1)) == 0_pInt) write(777,'(A)') ''
 enddo
 close(777)
 
!--------------------------------------------------------------------------------------------------
! set F to average values
 select case(myspectralsolver)
   case('basic')
     allocate(spectralF33New(3,3,resNew(1),resNew(2),resNew(3)))
     spectralF33New = spread(spread(spread(Favg,3,resNew(1)),4,resNew(2)),5,resNew(3))
     call IO_write_jobBinaryFile(777,'F',size(spectralF33New))
     write (777,rec=1) spectralF33New
     close (777)
     
   case('basicpetsc','al')
     allocate(spectralF9New(9,resNew(1),resNew(2),resNew(3)))
     spectralF9New = spread(spread(spread(reshape(Favg,[9]),2,resNew(1)),3,resNew(2)),4,resNew(3))
     call IO_write_jobBinaryFile(777,'F',size(spectralF9New))
     write (777,rec=1) spectralF9New
     close (777)
  end select

!---------------------------------------------------------------------------------
 allocate(F_lastIncNew(3,3,resNew(1),resNew(2),resNew(3)))

 call IO_read_jobBinaryFile(777,'F_aim_lastInc', &
                                 trim(getSolverJobName()),size(Favg_LastInc))
 read (777,rec=1) Favg_LastInc
 close (777)
 
 F_lastIncNew = spread(spread(spread(Favg_LastInc,3,resNew(1)),4,resNew(2)),5,resNew(3))
 
 call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',size(F_LastIncNew))
 write (777,rec=1) F_LastIncNew
 close (777)
 
 deallocate(F_lastIncNew)

! relocating data of material subroutine ---------------------------------------------------------
 allocate(material_phase    (1,1, mesh_NcpElems))
 allocate(material_phaseNew (1,1, NpointsNew))
 call IO_read_jobBinaryIntFile(777,'recordedPhase',trim(getSolverJobName()),size(material_phase))
 read (777,rec=1) material_phase
 close (777)
 do i = 1, NpointsNew
   material_phaseNew(1,1,i) = material_phase(1,1,indices(i))
 enddo
 do i = 1, mesh_NcpElems
   if (all(material_phaseNew(1,1,:) /= material_phase(1,1,i))) then
     write(6,*) 'mismatch in regridding'
     write(6,*) material_phase(1,1,i), 'not found in material_phaseNew'
   endif
 enddo
 call IO_write_jobBinaryIntFile(777,'recordedPhase',size(material_phaseNew))
 write (777,rec=1) material_phaseNew
 close (777) 
 deallocate(material_phase)
 deallocate(material_phaseNew)
!---------------------------------------------------------------------------
 allocate(F    (3,3,1,1, mesh_NcpElems))
 allocate(FNew (3,3,1,1, NpointsNew))
 call IO_read_jobBinaryFile(777,'convergedF',trim(getSolverJobName()),size(F))
 read (777,rec=1) F
 close (777)
 do i = 1, NpointsNew
   FNew(1:3,1:3,1,1,i) = F(1:3,1:3,1,1,indices(i))
 enddo

 call IO_write_jobBinaryFile(777,'convergedF',size(FNew))
 write (777,rec=1) FNew
 close (777) 
 deallocate(F)
 deallocate(FNew)
!--------------------------------------------------------------------- 
 allocate(Fp       (3,3,1,1,mesh_NcpElems))
 allocate(FpNew (3,3,1,1,NpointsNew))  
 call IO_read_jobBinaryFile(777,'convergedFp',trim(getSolverJobName()),size(Fp))
 read (777,rec=1) Fp
 close (777) 
 do i = 1, NpointsNew
   FpNew(1:3,1:3,1,1,i) = Fp(1:3,1:3,1,1,indices(i))
 enddo
 
 call IO_write_jobBinaryFile(777,'convergedFp',size(FpNew))
 write (777,rec=1) FpNew
 close (777) 
 deallocate(Fp)
 deallocate(FpNew)
!------------------------------------------------------------------------
 allocate(Lp       (3,3,1,1,mesh_NcpElems))
 allocate(LpNew  (3,3,1,1,NpointsNew)) 
 call IO_read_jobBinaryFile(777,'convergedLp',trim(getSolverJobName()),size(Lp))
 read (777,rec=1) Lp
 close (777)
 do i = 1, NpointsNew
   LpNew(1:3,1:3,1,1,i) = Lp(1:3,1:3,1,1,indices(i))
 enddo
 call IO_write_jobBinaryFile(777,'convergedLp',size(LpNew))
 write (777,rec=1) LpNew
 close (777)
 deallocate(Lp)
 deallocate(LpNew)
!----------------------------------------------------------------------------
 allocate(dcsdE       (6,6,1,1,mesh_NcpElems)) 
 allocate(dcsdENew (6,6,1,1,NpointsNew)) 
 call IO_read_jobBinaryFile(777,'convergeddcsdE',trim(getSolverJobName()),size(dcsdE))
 read (777,rec=1) dcsdE
 close (777)
 do i = 1, NpointsNew
   dcsdENew(1:6,1:6,1,1,i) = dcsdE(1:6,1:6,1,1,indices(i))
 enddo
 call IO_write_jobBinaryFile(777,'convergeddcsdE',size(dcsdENew))
 write (777,rec=1) dcsdENew
 close (777)
 deallocate(dcsdE)
 deallocate(dcsdENew)
!---------------------------------------------------------------------------
 allocate(dPdF       (3,3,3,3,1,1,mesh_NcpElems))
 allocate(dPdFNew (3,3,3,3,1,1,NpointsNew)) 
 call IO_read_jobBinaryFile(777,'convergeddPdF',trim(getSolverJobName()),size(dPdF))
 read (777,rec=1) dPdF
 close (777)
 do i = 1, NpointsNew
   dPdFNew(1:3,1:3,1:3,1:3,1,1,i) = dPdF(1:3,1:3,1:3,1:3,1,1,indices(i))
 enddo
 call IO_write_jobBinaryFile(777,'convergeddPdF',size(dPdFNew))
 write (777,rec=1) dPdFNew
 close (777)
 deallocate(dPdF)
 deallocate(dPdFNew)
!---------------------------------------------------------------------------
 allocate(Tstar        (6,1,1,mesh_NcpElems))
 allocate(TstarNew  (6,1,1,NpointsNew)) 
 call IO_read_jobBinaryFile(777,'convergedTstar',trim(getSolverJobName()),size(Tstar))
 read (777,rec=1) Tstar
 close (777)
 do i = 1, NpointsNew
   TstarNew(1:6,1,1,i) = Tstar(1:6,1,1,indices(i))
 enddo
 call IO_write_jobBinaryFile(777,'convergedTstar',size(TstarNew))
 write (777,rec=1) TstarNew
 close (777)
 deallocate(Tstar)
 deallocate(TstarNew)
 
! for the state, we first have to know the size------------------------------------------------------------------ 
 allocate(sizeStateConst(1,1,mesh_NcpElems))
 call IO_read_jobBinaryIntFile(777,'sizeStateConst',trim(getSolverJobName()),size(sizeStateConst))
 read (777,rec=1) sizeStateConst
 close (777)
 maxsize = maxval(sizeStateConst(1,1,1:mesh_NcpElems))
 allocate(StateConst      (1,1,mesh_NcpElems,maxsize))

 call IO_read_jobBinaryFile(777,'convergedStateConst',trim(getSolverJobName()))
 k = 0_pInt
 do i =1, mesh_NcpElems
   do j = 1,sizeStateConst(1,1,i)
     k = k+1_pInt
     read(777,rec=k) StateConst(1,1,i,j)
   enddo
 enddo
 close(777)
 call IO_write_jobBinaryFile(777,'convergedStateConst')
 k = 0_pInt
 do i = 1,NpointsNew
   do j = 1,sizeStateConst(1,1,indices(i))
     k=k+1_pInt
     write(777,rec=k) StateConst(1,1,indices(i),j)
   enddo
 enddo
 close (777)
 deallocate(sizeStateConst)
 deallocate(StateConst)
!---------------------------------------------------------------------------- 
 allocate(sizeStateHomog(1,mesh_NcpElems))
 call IO_read_jobBinaryIntFile(777,'sizeStateHomog',trim(getSolverJobName()),size(sizeStateHomog))
 read (777,rec=1) sizeStateHomog
 close (777)
 maxsize = maxval(sizeStateHomog(1,1:mesh_NcpElems))
 allocate(stateHomog      (1,mesh_NcpElems,maxsize))

 call IO_read_jobBinaryFile(777,'convergedStateHomog',trim(getSolverJobName()))
 k = 0_pInt
 do i =1, mesh_NcpElems
   do j = 1,sizeStateHomog(1,i)
     k = k+1_pInt
     read(777,rec=k) stateHomog(1,i,j)
   enddo
 enddo
 close(777)
 call IO_write_jobBinaryFile(777,'convergedStateHomog')
 k = 0_pInt
 do i = 1,NpointsNew
   do j = 1,sizeStateHomog(1,indices(i))
     k=k+1_pInt
     write(777,rec=k) stateHomog(1,indices(i),j)
   enddo
 enddo
 close (777)
 deallocate(sizeStateHomog)
 deallocate(stateHomog)
  
 deallocate(indices)
 write(6,*) 'finished regridding'
 
end function mesh_regrid


!--------------------------------------------------------------------------------------------------
!> @brief builds mesh of (distorted) cubes for given coordinates (= center of the cubes)
!--------------------------------------------------------------------------------------------------
function mesh_nodesAroundCentres(gDim,Favg,centres) result(nodes)
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic
 use math, only: &
   math_mul33x3
 
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
       if (k==0_pInt .or. k==iRes(3)+1_pInt .or. &                                                   ! z skin
           j==0_pInt .or. j==iRes(2)+1_pInt .or. &                                                   ! y skin
           i==0_pInt .or. i==iRes(1)+1_pInt      ) then                                              ! x skin
         me = [i,j,k]                                                                                ! me on skin
         shift = sign(abs(iRes+diag-2_pInt*me)/(iRes+diag),iRes+diag-2_pInt*me)
         lookup = me-diag+shift*iRes
         wrappedCentres(1:3,i+1_pInt,        j+1_pInt,        k+1_pInt) = &
                centres(1:3,lookup(1)+1_pInt,lookup(2)+1_pInt,lookup(3)+1_pInt) - &
                                                           math_mul33x3(Favg, shift*gDim)
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
!> @brief calculate coordinates in current configuration for given defgrad using linear 
! interpolation
!--------------------------------------------------------------------------------------------------
function mesh_deformedCoordsLinear(gDim,F,FavgIn) result(coords)
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic
 use math, only: &
   math_mul33x3

 implicit none
 real(pReal),   intent(in), dimension(:,:,:,:,:) :: &
   F
 real(pReal),               dimension(3,size(F,3),size(F,4),size(F,5)) :: &
   coords
 real(pReal),   intent(in), dimension(3) :: &
   gDim
 real(pReal),   intent(in), dimension(3,3), optional  :: &
   FavgIn
 real(pReal),               dimension(3,0:size(F,3)-1,0:size(F,4)-1,0:size(F,5)-1,0:7) :: &
   coordsAvgOrder
 integer(pInt), parameter,  dimension(3) ::  &
   iOnes = 1_pInt
 real(pReal),   parameter,  dimension(3) ::  &
   fOnes = 1.0_pReal
 real(pReal),               dimension(3) ::  &
   myStep, &
   negative, &
   positive, &
   offsetCoords, &
   parameterCoords, &
   stepLength, &
   fRes
 real(pReal),               dimension(3,3) :: &
   Favg
 integer(pInt), dimension(3) :: &
   rear, &
   init, &
   oppo, &
   me, &
   smallRes, &
   iRes
 integer(pInt) :: &
   i, j, k, s, o
 integer(pInt), parameter, dimension(3,0:7) :: &
   corner = reshape([ &
                                              0_pInt, 0_pInt, 0_pInt,&
                                              1_pInt, 0_pInt, 0_pInt,&
                                              1_pInt, 1_pInt, 0_pInt,&
                                              0_pInt, 1_pInt, 0_pInt,&
                                              1_pInt, 1_pInt, 1_pInt,&
                                              0_pInt, 1_pInt, 1_pInt,&
                                              0_pInt, 0_pInt, 1_pInt,&
                                              1_pInt, 0_pInt, 1_pInt &
                                                  ],[3,8]), &
   step = reshape([&
                                            1_pInt, 1_pInt, 1_pInt,&
                                           -1_pInt, 1_pInt, 1_pInt,&
                                           -1_pInt,-1_pInt, 1_pInt,&
                                            1_pInt,-1_pInt, 1_pInt,&
                                           -1_pInt,-1_pInt,-1_pInt,&
                                            1_pInt,-1_pInt,-1_pInt,&
                                            1_pInt, 1_pInt,-1_pInt,&
                                           -1_pInt, 1_pInt,-1_pInt &
                                                ], [3,8])
 integer(pInt), parameter, dimension(3,6) :: &
   order = reshape([ &
                                            1_pInt, 2_pInt, 3_pInt,&
                                            1_pInt, 3_pInt, 2_pInt,&
                                            2_pInt, 1_pInt, 3_pInt,&
                                            2_pInt, 3_pInt, 1_pInt,&
                                            3_pInt, 1_pInt, 2_pInt,&
                                            3_pInt, 2_pInt, 1_pInt &
                                                ], [3,6])
 
!--------------------------------------------------------------------------------------------------
! initializing variables
 iRes      =  [size(F,3),size(F,4),size(F,5)]
 fRes      =  real(iRes,pReal)
 smallRes =  iRes - 1_pInt
 coordsAvgOrder = 0.0_pReal
 stepLength = gDim/fRes

!--------------------------------------------------------------------------------------------------
! report
 if (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt) then
   write(6,'(a)')          ' Restore geometry using linear integration'
   write(6,'(a,3(i12  ))') ' grid     a b c: ', iRes
   write(6,'(a,3(f12.5))') ' size     x y z: ', gDim
 endif

!--------------------------------------------------------------------------------------------------
! determine average deformation gradient
 if (present(FavgIn)) then
   if (all(FavgIn < 0.0_pReal)) then                                                                ! the f2py way to tell it is not present
     Favg = sum(sum(sum(F,dim=5),dim=4),dim=3) / product(fRes)
       else
     Favg = FavgIn
   endif
 else
   Favg = sum(sum(sum(F,dim=5),dim=4),dim=3) / product(fRes)
 endif

!--------------------------------------------------------------------------------------------------
! loop over starting corners (from 0 to 7)
 cornerLooping: do s = 0_pInt, 7_pInt
   init = corner(1:3,s) * smallRes  + iOnes
   oppo = corner(1:3,mod((s+4_pInt),8_pInt)) * smallRes + iOnes

!--------------------------------------------------------------------------------------------------
! permutation of ways on each corner
   permutationLooping: do o = 1_pInt,6_pInt
     coords = 0.0_pReal
     do k = init(order(3,o)), oppo(order(3,o)), step(order(3,o),s)
       rear(order(2,o)) = init(order(2,o))
       do j = init(order(2,o)), oppo(order(2,o)), step(order(2,o),s)
         rear(order(1,o)) = init(order(1,o))
         do i = init(order(1,o)), oppo(order(1,o)), step(order(1,o),s)
           me(order(1:3,o)) = [i,j,k]
           if ( all(me == init)) then
             coords(1:3,me(1),me(2),me(3)) = gDim*( math_mul33x3(Favg,real(corner(1:3,s),pReal)) &
                                                      +math_mul33x3(F(1:3,1:3,me(1),me(2),me(3)), &
                                                          0.5_pReal*real(step(1:3,s)/iRes,pReal)))
           else
             myStep = (me-rear)*stepLength
             coords(1:3,me(1),me(2),me(3)) = coords(1:3,rear(1),rear(2),rear(3)) + &
                                           0.5_pReal*math_mul33x3(F(1:3,1:3,me(1),me(2),me(3)) + &
                                                   F(1:3,1:3,rear(1),rear(2),rear(3)),myStep)
           endif
           rear = me
     enddo; enddo; enddo
     coordsAvgOrder(1:3,0:smallRes(1),0:smallRes(2),0:smallRes(3),s) = &
         coordsAvgOrder(1:3,0:smallRes(1),0:smallRes(2),0:smallRes(3),s) + coords/6.0_pReal
   enddo permutationLooping
   offsetCoords = coordsAvgOrder(1:3,0,0,0,s)
   do k = 0_pInt, smallRes(3); do j = 0_pInt, smallRes(2); do i = 0_pInt, smallRes(1)
     coordsAvgOrder(1:3,i,j,k,s) = coordsAvgOrder(1:3,i,j,k,s) - offsetCoords
   enddo; enddo; enddo
 enddo cornerLooping

!--------------------------------------------------------------------------------------------------
! linear interpolation starting at each corner (comparable to linear shape function FEM)
 do k = 0_pInt, smallRes(3); do j = 0_pInt, smallRes(2); do i = 0_pInt, smallRes(1)
   parameterCoords = (2.0_pReal*real([i,j,k]+1,pReal)-fRes)/fRes
   positive = fones + parameterCoords
   negative = fones - parameterCoords
   coords(1:3,i+1_pInt,j+1_pInt,k+1_pInt) &
              =(coordsAvgOrder(1:3,i,j,k,0) *negative(1)*negative(2)*negative(3)&
              + coordsAvgOrder(1:3,i,j,k,1) *positive(1)*negative(2)*negative(3)&
              + coordsAvgOrder(1:3,i,j,k,2) *positive(1)*positive(2)*negative(3)&
              + coordsAvgOrder(1:3,i,j,k,3) *negative(1)*positive(2)*negative(3)&
              + coordsAvgOrder(1:3,i,j,k,4) *positive(1)*positive(2)*positive(3)&
              + coordsAvgOrder(1:3,i,j,k,5) *negative(1)*positive(2)*positive(3)&
              + coordsAvgOrder(1:3,i,j,k,6) *negative(1)*negative(2)*positive(3)&
              + coordsAvgOrder(1:3,i,j,k,7) *positive(1)*negative(2)*positive(3))*0.125_pReal
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! setting base node to (0,0,0)
 offsetCoords = math_mul33x3(F(1:3,1:3,1,1,1),stepLength/2.0_pReal) - coords(1:3,1,1,1)
 do k = 1_pInt, iRes(3); do j = 1_pInt, iRes(2); do i = 1_pInt, iRes(1)
   coords(1:3,i,j,k) = coords(1:3,i,j,k) + offsetCoords
 enddo; enddo; enddo

end function mesh_deformedCoordsLinear


!--------------------------------------------------------------------------------------------------
!> @brief calculate coordinates in current configuration for given defgrad
! using integration in Fourier space 
!--------------------------------------------------------------------------------------------------
function mesh_deformedCoordsFFT(gDim,F,FavgIn,scalingIn) result(coords)
 use IO, only: &
   IO_error
 use numerics, only: &
   fftw_timelimit, &
   fftw_planner_flag
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic
 use math, only: &
   PI, &
   math_mul33x3

 implicit none
 real(pReal), intent(in), dimension(:,:,:,:,:) ::                                F
 real(pReal),             dimension(3,size(F,3),size(F,4),size(F,5)) ::          coords
 real(pReal), intent(in), dimension(3) ::                                        gDim
 real(pReal), intent(in), dimension(3,3),                           optional  :: FavgIn
 real(pReal), intent(in), dimension(3),                             optional  :: scalingIn

! allocatable arrays for fftw c routines
 type(C_PTR) :: planForth, planBack
 type(C_PTR) :: coords_fftw, defgrad_fftw
 real(pReal),    dimension(:,:,:,:,:), pointer :: F_real
 complex(pReal), dimension(:,:,:,:,:), pointer :: F_fourier
 real(pReal),    dimension(:,:,:,:),   pointer :: coords_real
 complex(pReal), dimension(:,:,:,:),   pointer :: coords_fourier
 ! other variables
 integer(pInt) :: i, j, k, m, res1Red
 integer(pInt), dimension(3) :: k_s, iRes
 real(pReal),   dimension(3) :: scaling, step, offset_coords, integrator
 real(pReal),   dimension(3,3) :: Favg
 integer(pInt), dimension(2:3,2) :: Nyquist                                                         ! highest frequencies to be removed (1 if even, 2 if odd)

 if (present(scalingIn)) then
   where (scalingIn < 0.0_pReal)                                                                    ! the f2py way to tell it is not present
     scaling = [1.0_pReal,1.0_pReal,1.0_pReal]
   elsewhere
     scaling = scalingIn
   endwhere
 else
   scaling = 1.0_pReal
 endif
 
 iRes =  [size(F,3),size(F,4),size(F,5)]
 integrator = gDim / 2.0_pReal / PI                                                                 ! see notes where it is used
 res1Red = iRes(1)/2_pInt + 1_pInt                                                                 ! size of complex array in first dimension (c2r, r2c)
 step = gDim/real(iRes, pReal)
 Nyquist(2,1:2) = [res(2)/2_pInt + 1_pInt, res(2)/2_pInt + 1_pInt + mod(res(2),2_pInt)]
 Nyquist(3,1:2) = [res(3)/2_pInt + 1_pInt, res(3)/2_pInt + 1_pInt + mod(res(3),2_pInt)]

!--------------------------------------------------------------------------------------------------
! report
 if (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt) then
   write(6,'(a)')          ' Restore geometry using FFT-based integration'
   write(6,'(a,3(i12  ))') ' grid     a b c: ', iRes
   write(6,'(a,3(f12.5))') ' size     x y z: ', gDim
 endif

!--------------------------------------------------------------------------------------------------
! sanity check
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) &
   call IO_error(0_pInt,ext_msg='Fortran to C in mesh_deformedCoordsFFT')

!--------------------------------------------------------------------------------------------------
! allocation and FFTW initialization
 defgrad_fftw =  fftw_alloc_complex(int(res1Red *iRes(2)*iRes(3)*9_pInt,C_SIZE_T))                 ! C_SIZE_T is of type integer(8)
 coords_fftw  =  fftw_alloc_complex(int(res1Red *iRes(2)*iRes(3)*3_pInt,C_SIZE_T))                 ! C_SIZE_T is of type integer(8)
 call c_f_pointer(defgrad_fftw, F_real, &
   [iRes(1)+2_pInt-mod(iRes(1),2_pInt),iRes(2),iRes(3),3_pInt,3_pInt])
 call c_f_pointer(defgrad_fftw, F_fourier, &
   [res1Red,                           iRes(2),iRes(3),3_pInt,3_pInt])
 call c_f_pointer(coords_fftw, coords_real, &
   [iRes(1)+2_pInt-mod(iRes(1),2_pInt),iRes(2),iRes(3),3_pInt])
 call c_f_pointer(coords_fftw, coords_fourier, &
   [res1Red,                           iRes(2),iRes(3),3_pInt])

 call fftw_set_timelimit(fftw_timelimit)
 planForth = fftw_plan_many_dft_r2c(3_pInt,[iRes(3),iRes(2) ,iRes(1)],9_pInt,&                      ! dimensions , length in each dimension in reversed order
                                    F_real,[iRes(3),iRes(2) ,iRes(1)+2_pInt-mod(iRes(1),2_pInt)],&  ! input data , physical length in each dimension in reversed order
                                    1_pInt, iRes(3)*iRes(2)*(iRes(1)+2_pInt-mod(iRes(1),2_pInt)),&  ! striding   , product of physical lenght in the 3 dimensions
                                 F_fourier,[iRes(3),iRes(2) ,res1Red],&
                                   1_pInt,  iRes(3)*iRes(2)* res1Red,fftw_planner_flag)

 planBack  = fftw_plan_many_dft_c2r(3_pInt,[iRes(3),iRes(2) ,iRes(1)],3_pInt,&
                            coords_fourier,[iRes(3),iRes(2) ,res1Red],&
                                   1_pInt,  iRes(3)*iRes(2)* res1Red,&
                               coords_real,[iRes(3),iRes(2) ,iRes(1)+2_pInt-mod(iRes(1),2_pInt)],&
                                   1_pInt,  iRes(3)*iRes(2)*(iRes(1)+2_pInt-mod(iRes(1),2_pInt)),&
                                                                      fftw_planner_flag)
 F_real(1:iRes(1),1:iRes(2),1:iRes(3),1:3,1:3) = & 
                                         reshape(F,[res(1),res(2),res(3),3,3], order = [4,5,1,2,3]) ! F_real is overwritten during plan creatio, is larger (padding) and has different order

!--------------------------------------------------------------------------------------------------
! FFT
 call fftw_execute_dft_r2c(planForth, F_real, F_fourier)
 
!--------------------------------------------------------------------------------------------------
! if no average F is given, compute it in Fourier space
 if (present(FavgIn)) then
   if (all(FavgIn < 0.0_pReal)) then
     Favg = real(F_fourier(1,1,1,1:3,1:3),pReal)/real(product(iRes),pReal)                          !the f2py way to tell it is not present
   else
     Favg = FavgIn
   endif
 else
   Favg = real(F_fourier(1,1,1,1:3,1:3),pReal)/real(product(iRes),pReal)
 endif
 
!--------------------------------------------------------------------------------------------------
! remove highest frequency in each direction, in third direction only if not 2D

 if(iRes(1)/=1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
   F_fourier (res1Red,  1:iRes(2),                1:iRes(3),              1:3,1:3) &
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
 if(iRes(2)/=1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
   F_fourier (1:res1Red,Nyquist(2,1):Nyquist(2,2),1:iRes(3),               1:3,1:3) & 
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
 if(iRes(3)/=1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
   F_fourier (1:res1Red,1:iRes(2),                Nyquist(3,1):Nyquist(3,2),1:3,1:3) &
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)

!--------------------------------------------------------------------------------------------------
! integration in Fourier space
 coords_fourier = cmplx(0.0_pReal,0.0_pReal,pReal)
 do k = 1_pInt, iRes(3)
   k_s(3) = k-1_pInt
   if(k > iRes(3)/2_pInt+1_pInt) k_s(3) = k_s(3)-iRes(3)
   do j = 1_pInt, iRes(2)
     k_s(2) = j-1_pInt
     if(j > iRes(2)/2_pInt+1_pInt) k_s(2) = k_s(2)-iRes(2)
     do i = 1_pInt, res1Red
       k_s(1) = i-1_pInt
       do m = 1_pInt,3_pInt
         coords_fourier(i,j,k,m) = sum(F_fourier(i,j,k,m,1:3)*&
                                       cmplx(0.0_pReal,real(k_s,pReal)*integrator,pReal))
       enddo
       if (any(k_s /= 0_pInt)) coords_fourier(i,j,k,1:3) = &
                               coords_fourier(i,j,k,1:3) / cmplx(-sum(k_s*k_s),0.0_pReal,pReal)
 enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! iFFT and freeing memory
 call fftw_execute_dft_c2r(planBack,coords_fourier,coords_real)
 coords = reshape(coords_real(1:iRes(1),1:iRes(2),1:iRes(3),1:3), [3,iRes(1),iRes(2),iRes(3)], &
                  order = [2,3,4,1])/real(product(iRes),pReal)                                      ! weight and change order
 call fftw_destroy_plan(planForth)
 call fftw_destroy_plan(planBack)
 call fftw_free(defgrad_fftw)
 call fftw_free(coords_fftw)

!--------------------------------------------------------------------------------------------------
! add average to scaled fluctuation and put (0,0,0) on (0,0,0)
 offset_coords = math_mul33x3(F(1:3,1:3,1,1,1),step/2.0_pReal) - scaling*coords(1:3,1,1,1)
 forall(k = 1_pInt:iRes(3), j = 1_pInt:iRes(2), i = 1_pInt:iRes(1)) &
   coords(1:3,i,j,k) = scaling(1:3)*coords(1:3,i,j,k) &
                     + offset_coords &
                     + math_mul33x3(Favg,step*real([i,j,k]-1_pInt,pReal))

end function mesh_deformedCoordsFFT


!--------------------------------------------------------------------------------------------------
!> @brief calculates the mismatch between volume of reconstructed (compatible) cube and 
! determinant of defgrad at the FP
!--------------------------------------------------------------------------------------------------
function mesh_volumeMismatch(gDim,F,nodes) result(vMismatch)
 use IO, only: &
   IO_error
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic
 use math, only:  &
   PI, &
   math_det33, &
   math_volTetrahedron

 implicit none
 real(pReal),   intent(in), dimension(:,:,:,:,:) :: &
   F
 real(pReal),               dimension(size(F,3),size(F,4),size(F,5)) :: &
   vMismatch
 real(pReal),   intent(in), dimension(:,:,:,:)   :: &
   nodes
 real(pReal),               dimension(3) :: &
   gDim
 integer(pInt),             dimension(3) :: &
   iRes
 real(pReal),   dimension(3,8) ::  coords
 integer(pInt) :: i,j,k
 real(pReal) :: volInitial

 iRes = [size(F,3),size(F,4),size(F,5)]
 volInitial = product(gDim)/real(product(iRes), pReal)

!--------------------------------------------------------------------------------------------------
! report and check 
 if (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt) then
   write(6,'(a)')          ' Calculating volume mismatch'
   write(6,'(a,3(i12  ))') ' grid     a b c: ', iRes
   write(6,'(a,3(f12.5))') ' size     x y z: ', gDim
 endif
  
 if (any([iRes/=size(nodes,2)-1_pInt,iRes/=size(nodes,3)-1_pInt,iRes/=size(nodes,4)-1_pInt]))&
   call IO_error(0_pInt,ext_msg='Arrays F and nodes in mesh_volumeMismatch')
 
!--------------------------------------------------------------------------------------------------
! calculate actual volume and volume resulting from deformation gradient
 do k = 1_pInt,iRes(3)
   do j = 1_pInt,iRes(2)
     do i = 1_pInt,iRes(1)
       coords(1:3,1) = nodes(1:3,i,       j,       k       )
       coords(1:3,2) = nodes(1:3,i+1_pInt,j,       k       )
       coords(1:3,3) = nodes(1:3,i+1_pInt,j+1_pInt,k       )
       coords(1:3,4) = nodes(1:3,i,       j+1_pInt,k       )
       coords(1:3,5) = nodes(1:3,i,       j,       k+1_pInt)
       coords(1:3,6) = nodes(1:3,i+1_pInt,j,       k+1_pInt)
       coords(1:3,7) = nodes(1:3,i+1_pInt,j+1_pInt,k+1_pInt)
       coords(1:3,8) = nodes(1:3,i,       j+1_pInt,k+1_pInt)
       vMismatch(i,j,k) = &
           abs(math_volTetrahedron(coords(1:3,7),coords(1:3,1),coords(1:3,8),coords(1:3,4))) &
         + abs(math_volTetrahedron(coords(1:3,7),coords(1:3,1),coords(1:3,8),coords(1:3,5))) &
         + abs(math_volTetrahedron(coords(1:3,7),coords(1:3,1),coords(1:3,3),coords(1:3,4))) &
         + abs(math_volTetrahedron(coords(1:3,7),coords(1:3,1),coords(1:3,3),coords(1:3,2))) &
         + abs(math_volTetrahedron(coords(1:3,7),coords(1:3,5),coords(1:3,2),coords(1:3,6))) &
         + abs(math_volTetrahedron(coords(1:3,7),coords(1:3,5),coords(1:3,2),coords(1:3,1)))
       vMismatch(i,j,k) = vMismatch(i,j,k)/math_det33(F(1:3,1:3,i,j,k))
 enddo; enddo; enddo
 vMismatch = vMismatch/volInitial

end function mesh_volumeMismatch


!--------------------------------------------------------------------------------------------------
!> @brief Routine to calculate the mismatch between the vectors from the central point to
! the corners of reconstructed (combatible) volume element and the vectors calculated by deforming
! the initial volume element with the  current deformation gradient
!--------------------------------------------------------------------------------------------------
function mesh_shapeMismatch(gDim,F,nodes,centres)  result(sMismatch)
 use IO, only: &
   IO_error
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic
 use math, only: &
  math_mul33x3

 implicit none
 real(pReal),   intent(in), dimension(:,:,:,:,:) :: &
   F
 real(pReal),               dimension(size(F,3),size(F,4),size(F,5)) :: &
   sMismatch
 real(pReal),   intent(in), dimension(:,:,:,:)   :: &
   nodes, &
   centres
 real(pReal),               dimension(3) :: &
   gDim, &
   fRes
 integer(pInt),             dimension(3) :: &
   iRes
 real(pReal), dimension(3,8) :: coordsInitial
 integer(pInt) i,j,k

 iRes = [size(F,3),size(F,4),size(F,5)]
 fRes = real(iRes,pReal)
 
!--------------------------------------------------------------------------------------------------
! report and check
 if (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt) then
   write(6,'(a)')          ' Calculating shape mismatch'
   write(6,'(a,3(i12  ))') ' grid     a b c: ', iRes
   write(6,'(a,3(f12.5))') ' size     x y z: ', gDim
 endif

 if(any([iRes/=size(nodes,2)-1_pInt,iRes/=size(nodes,3)-1_pInt,iRes/=size(nodes,4)-1_pInt]) .or.&
    any([iRes/=size(centres,2),     iRes/=size(centres,3),     iRes/=size(centres,4)]))&
   call IO_error(0_pInt,ext_msg='Arrays F and nodes/centres in mesh_shapeMismatch')
   
!--------------------------------------------------------------------------------------------------
! initial positions
 coordsInitial(1:3,1) = [-gDim(1)/fRes(1),-gDim(2)/fRes(2),-gDim(3)/fRes(3)]
 coordsInitial(1:3,2) = [+gDim(1)/fRes(1),-gDim(2)/fRes(2),-gDim(3)/fRes(3)]
 coordsInitial(1:3,3) = [+gDim(1)/fRes(1),+gDim(2)/fRes(2),-gDim(3)/fRes(3)]
 coordsInitial(1:3,4) = [-gDim(1)/fRes(1),+gDim(2)/fRes(2),-gDim(3)/fRes(3)]
 coordsInitial(1:3,5) = [-gDim(1)/fRes(1),-gDim(2)/fRes(2),+gDim(3)/fRes(3)]
 coordsInitial(1:3,6) = [+gDim(1)/fRes(1),-gDim(2)/fRes(2),+gDim(3)/fRes(3)]
 coordsInitial(1:3,7) = [+gDim(1)/fRes(1),+gDim(2)/fRes(2),+gDim(3)/fRes(3)]
 coordsInitial(1:3,8) = [-gDim(1)/fRes(1),+gDim(2)/fRes(2),+gDim(3)/fRes(3)]
 coordsInitial = coordsInitial/2.0_pReal
 
!--------------------------------------------------------------------------------------------------
! compare deformed original and deformed positions to actual positions
 do k = 1_pInt,iRes(3)
   do j = 1_pInt,iRes(2)
     do i = 1_pInt,iRes(1)
       sMismatch(i,j,k) = &
           sqrt(sum((nodes(1:3,i,       j,       k        ) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,1)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i+1_pInt,j,       k        ) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,2)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i+1_pInt,j+1_pInt,k        ) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,3)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i,       j+1_pInt,k        ) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,4)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i,       j,       k+1_pInt) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,5)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i+1_pInt,j,       k+1_pInt) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,6)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i+1_pInt,j+1_pInt,k+1_pInt) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,7)))**2.0_pReal))&
         + sqrt(sum((nodes(1:3,i,       j+1_pInt,k+1_pInt) - centres(1:3,i,j,k)&
                    - math_mul33x3(F(1:3,1:3,i,j,k), coordsInitial(1:3,8)))**2.0_pReal))
 enddo; enddo; enddo

end function mesh_shapeMismatch
#endif


#ifdef Marc
!--------------------------------------------------------------------------------------------------
!> @brief Figures out table styles (Marc only) and stores to 'initialcondTableStyle' and 
!! 'hypoelasticTableStyle'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_get_tableStyles(myUnit)
 use IO, only: &
   IO_lc, &
   IO_intValue, &
   IO_stringValue, &
   IO_stringPos
 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 6_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line

 initialcondTableStyle = 0_pInt
 hypoelasticTableStyle = 0_pInt
 
610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'table' .and. myPos(1_pInt) > 5) then
     initialcondTableStyle = IO_intValue(line,myPos,4_pInt)
     hypoelasticTableStyle = IO_intValue(line,myPos,5_pInt)
     exit
   endif
 enddo

620 end subroutine mesh_marc_get_tableStyles


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of nodes and elements in mesh and stores them in
!! 'mesh_Nelems' and 'mesh_Nnodes'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_count_nodesAndElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_IntValue
 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,myPos,1_pInt)) == 'sizing') &
       mesh_Nelems = IO_IntValue (line,myPos,3_pInt)
   if ( IO_lc(IO_StringValue(line,myPos,1_pInt)) == 'coordinates') then
     read (myUnit,610,END=620) line
     myPos = IO_stringPos(line,maxNchunks)
     mesh_Nnodes = IO_IntValue (line,myPos,2_pInt)
     exit                                                                                          ! assumes that "coordinates" comes later in file
   endif
 enddo

620 end subroutine mesh_marc_count_nodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of element sets in mesh. Stores to 'mesh_NelemSets', and
!! 'mesh_maxNelemInSet'
!--------------------------------------------------------------------------------------------------
 subroutine mesh_marc_count_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countContinuousIntValues
                 
 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line

 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,myPos,1_pInt)) == 'define' .and. &
        IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'element' ) then
     mesh_NelemSets = mesh_NelemSets + 1_pInt
     mesh_maxNelemInSet = max(mesh_maxNelemInSet, &
                              IO_countContinuousIntValues(myUnit))
   endif
 enddo

620 end subroutine mesh_marc_count_elementSets


!********************************************************************
! map element sets
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!********************************************************************
subroutine mesh_marc_map_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_continuousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: elemSet = 0_pInt

 allocate (mesh_nameElemSet(mesh_NelemSets))                     ; mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1_pInt+mesh_maxNelemInSet,mesh_NelemSets)) ; mesh_mapElemSet = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=640) line
   myPos = IO_stringPos(line,maxNchunks)
   if( (IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'define' ) .and. &
       (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'element' ) ) then
      elemSet = elemSet+1_pInt
      mesh_nameElemSet(elemSet) = trim(IO_stringValue(line,myPos,4_pInt))
      mesh_mapElemSet(:,elemSet) = IO_continuousIntValues(myUnit,mesh_maxNelemInSet,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
   endif
 enddo
 
640 end subroutine mesh_marc_map_elementSets


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of CP elements in mesh and stores them in 'mesh_NcpElems'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_count_cpElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countContinuousIntValues
                 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 integer(pInt) :: i
 character(len=300):: line

 mesh_NcpElems = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'hypoelastic') then
       do i=1_pInt,3_pInt+hypoelasticTableStyle  ! Skip 3 or 4 lines
         read (myUnit,610,END=620) line
       enddo
       mesh_NcpElems = mesh_NcpElems + IO_countContinuousIntValues(myUnit)
     exit
   endif
 enddo

620 end subroutine mesh_marc_count_cpElements


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPelem'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_elements(myUnit)

 use math, only: math_qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_continuousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt), dimension (1_pInt+mesh_NcpElems) :: contInts
 integer(pInt) :: i,cpElem = 0_pInt

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=660) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'hypoelastic' ) then
     do i=1_pInt,3_pInt+hypoelasticTableStyle                                                       ! skip three (or four if new table style!) lines
       read (myUnit,610,END=660) line 
     enddo
     contInts = IO_continuousIntValues(myUnit,mesh_NcpElems,mesh_nameElemSet,&
                                              mesh_mapElemSet,mesh_NelemSets)
     do i = 1_pInt,contInts(1)
       cpElem = cpElem+1_pInt
       mesh_mapFEtoCPelem(1,cpElem) = contInts(1_pInt+i)
       mesh_mapFEtoCPelem(2,cpElem) = cpElem
     enddo
   endif
 enddo

660 call math_qsort(mesh_mapFEtoCPelem,1_pInt,int(size(mesh_mapFEtoCPelem,2_pInt),pInt))            ! should be mesh_NcpElems

end subroutine mesh_marc_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPnode'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_nodes(myUnit)

 use math, only: math_qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_fixedIntValue
                 
 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt) :: i

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

610 FORMAT(A300)

 node_count = 0_pInt

 rewind(myUnit)
 do
   read (myUnit,610,END=650) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'coordinates' ) then
     read (myUnit,610,END=650) line                                                                 ! skip crap line
     do i = 1_pInt,mesh_Nnodes
       read (myUnit,610,END=650) line
       mesh_mapFEtoCPnode(1_pInt,i) = IO_fixedIntValue (line,[ 0_pInt,10_pInt],1_pInt)
       mesh_mapFEtoCPnode(2_pInt,i) = i
     enddo
     exit
   endif
 enddo

650 call math_qsort(mesh_mapFEtoCPnode,1_pInt,int(size(mesh_mapFEtoCPnode,2_pInt),pInt))
 
end subroutine mesh_marc_map_nodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!! Allocates global arrays 'mesh_node0' and 'mesh_node'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_build_nodes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_fixedIntValue, &
                 IO_fixedNoEFloatValue
use numerics, only: numerics_unitlength

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), dimension(5), parameter :: node_ends = int([0,10,30,50,70],pInt)
 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,j,m

 allocate ( mesh_node0 (3,mesh_Nnodes) ); mesh_node0 = 0.0_pReal
 allocate ( mesh_node  (3,mesh_Nnodes) ); mesh_node  = 0.0_pReal

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=670) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'coordinates' ) then
     read (myUnit,610,END=670) line                                                                 ! skip crap line
     do i=1_pInt,mesh_Nnodes
       read (myUnit,610,END=670) line
       m = mesh_FEasCP('node',IO_fixedIntValue(line,node_ends,1_pInt))
       do j = 1_pInt,3_pInt
         mesh_node0(j,m) = numerics_unitlength * IO_fixedNoEFloatValue(line,node_ends,j+1_pInt)
       enddo
     enddo
     exit
   endif
 enddo

670 mesh_node = mesh_node0

end subroutine mesh_marc_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets maximum count of nodes, IPs, IP neighbors, and subNodes among cpElements.
!! Allocates global arrays 'mesh_maxNnodes', 'mesh_maxNips', mesh_maxNipNeighbors', 
!! and mesh_maxNsubNodes
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_count_cpSizes(myUnit)
 
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_skipChunks

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,t,g,e

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNsubNodes    = 0_pInt
 
610 FORMAT(A300)
 rewind(myUnit)
 do
   read (myUnit,610,END=630) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'connectivity' ) then
     read (myUnit,610,END=630) line                                                                 ! Garbage line
     do i=1_pInt,mesh_Nelems                                                                        ! read all elements
       read (myUnit,610,END=630) line
       myPos = IO_stringPos(line,maxNchunks)                                                        ! limit to id and type
       e = mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt))
       if (e /= 0_pInt) then
         t = FE_mapElemtype(IO_stringValue(line,myPos,2_pInt))
         g = FE_geomtype(t)
         mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(g))
         mesh_maxNips =         max(mesh_maxNips,FE_Nips(g))
         mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(g))
         mesh_maxNsubNodes =    max(mesh_maxNsubNodes,FE_NsubNodes(g))
         call IO_skipChunks(myUnit,FE_NoriginalNodes(t)-(myPos(1_pInt)-2_pInt))                     ! read on if FE_Nnodes exceeds node count present on current line
       endif
     enddo
     exit
   endif
 enddo
 
630 end subroutine mesh_marc_count_cpSizes


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, mat, tex, and node list per elemen.
!! Allocates global array 'mesh_element'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_build_elements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_fixedNoEFloatValue, &
                 IO_skipChunks, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_continuousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 66_pInt                                                  ! limit to 64 nodes max (plus ID, type)
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt), dimension(1_pInt+mesh_NcpElems) :: contInts
 integer(pInt) :: i,j,t,sv,myVal,e

 allocate (mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=620) line
   myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'connectivity' ) then
     read (myUnit,610,END=620) line                                                                 ! garbage line
     do i = 1_pInt,mesh_Nelems
       read (myUnit,610,END=620) line
       myPos = IO_stringPos(line,maxNchunks)
       e = mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt))
       if (e /= 0_pInt) then                                                                        ! disregard non CP elems
         t = FE_mapElemtype(IO_StringValue(line,myPos,2_pInt))                                      ! elem type
         mesh_element(2,e) = t
         mesh_element(1,e) = IO_IntValue (line,myPos,1_pInt)                                        ! FE id
           do j = 1_pInt,FE_Nnodes(FE_geomtype(t))
             mesh_element(j+4_pInt,e) = IO_IntValue(line,myPos,j+2_pInt)                            ! copy FE ids of nodes
           enddo  
           call IO_skipChunks(myUnit,FE_NoriginalNodes(t)-(myPos(1_pInt)-2_pInt))                   ! read on if FE_Nnodes exceeds node count present on current line
       endif
     enddo
     exit
   endif
 enddo
 
620 rewind(myUnit)                                                                                  ! just in case "initial state" apears before "connectivity"
 read (myUnit,610,END=620) line
 do
   myPos(1:1+2*2) = IO_stringPos(line,2_pInt)
   if( (IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'initial') .and. &
       (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'state') ) then
     if (initialcondTableStyle == 2_pInt) read (myUnit,610,END=620) line                            ! read extra line for new style
     read (myUnit,610,END=630) line                                                                 ! read line with index of state var
     myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
     sv = IO_IntValue(line,myPos,1_pInt)                                                            ! figure state variable index
     if( (sv == 2_pInt).or.(sv == 3_pInt) ) then                                                    ! only state vars 2 and 3 of interest
       read (myUnit,610,END=620) line                                                               ! read line with value of state var
       myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
       do while (scan(IO_stringValue(line,myPos,1_pInt),'+-',back=.true.)>1)                        ! is noEfloat value?
         myVal = nint(IO_fixedNoEFloatValue(line,[0_pInt,20_pInt],1_pInt),pInt)                     ! state var's value
         mesh_maxValStateVar(sv-1_pInt) = max(myVal,mesh_maxValStateVar(sv-1_pInt))                 ! remember max val of homogenization and microstructure index
         if (initialcondTableStyle == 2_pInt) then
           read (myUnit,610,END=630) line                                                           ! read extra line
           read (myUnit,610,END=630) line                                                           ! read extra line
         endif
         contInts = IO_continuousIntValues&                                                          ! get affected elements
                   (myUnit,mesh_Nelems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
         do i = 1_pInt,contInts(1)
           e = mesh_FEasCP('elem',contInts(1_pInt+i))
           mesh_element(1_pInt+sv,e) = myVal
         enddo
         if (initialcondTableStyle == 0_pInt) read (myUnit,610,END=620) line                        ! ignore IP range for old table style
         read (myUnit,610,END=630) line
         myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
       enddo
     endif
   else   
     read (myUnit,610,END=630) line
   endif
 enddo

630 end subroutine mesh_marc_build_elements
#endif 

#ifdef Abaqus
!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of nodes and elements in mesh and stores them in
!! 'mesh_Nelems' and 'mesh_Nnodes'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_nodesAndElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_error
                 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) :: line
 logical :: inPart

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt
 
610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.
   
   if (inPart .or. noPart) then
     select case ( IO_lc(IO_stringValue(line,myPos,1_pInt)))
       case('*node')
          if( &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'print'    .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'file'     .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' &
             ) &
            mesh_Nnodes = mesh_Nnodes + IO_countDataLines(myUnit)
       case('*element')
          if( &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'matrix'   .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' &
             ) then
            mesh_Nelems = mesh_Nelems + IO_countDataLines(myUnit)
          endif
     endselect
   endif
 enddo
 
620 if (mesh_Nnodes < 2_pInt)  call IO_error(error_ID=900_pInt)
 if (mesh_Nelems == 0_pInt) call IO_error(error_ID=901_pInt)
 
end subroutine mesh_abaqus_count_nodesAndElements


!********************************************************************
! count overall number of element sets in mesh
!
! mesh_NelemSets, mesh_maxNelemInSet
!********************************************************************
subroutine mesh_abaqus_count_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) :: line
 logical :: inPart
 
 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = mesh_Nelems               ! have to be conservative, since Abaqus allows for recursive definitons
 
610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.
   
   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*elset' ) &
     mesh_NelemSets = mesh_NelemSets + 1_pInt
 enddo

620 continue
 if (mesh_NelemSets == 0) call IO_error(error_ID=902_pInt)

end subroutine mesh_abaqus_count_elementSets


!********************************************************************
! count overall number of solid sections sets in mesh (Abaqus only)
!
! mesh_Nmaterials
!********************************************************************
subroutine mesh_abaqus_count_materials(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 logical inPart
 
 mesh_Nmaterials = 0_pInt
 
610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,myPos,1_pInt)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'section' ) &
     mesh_Nmaterials = mesh_Nmaterials + 1_pInt
 enddo

620 if (mesh_Nmaterials == 0_pInt) call IO_error(error_ID=903_pInt)
 
end subroutine mesh_abaqus_count_materials


!********************************************************************
! Build element set mapping 
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!********************************************************************
subroutine mesh_abaqus_map_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_continuousIntValues, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: elemSet = 0_pInt,i
 logical :: inPart = .false.

 allocate (mesh_nameElemSet(mesh_NelemSets))                          ; mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1_pInt+mesh_maxNelemInSet,mesh_NelemSets)) ; mesh_mapElemSet  = 0_pInt

610 FORMAT(A300)


 rewind(myUnit)
 do
   read (myUnit,610,END=640) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.
   
   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*elset' ) then
     elemSet = elemSet + 1_pInt
     mesh_nameElemSet(elemSet)  = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'elset'))
     mesh_mapElemSet(:,elemSet) = IO_continuousIntValues(myUnit,mesh_Nelems,mesh_nameElemSet,&
                                          mesh_mapElemSet,elemSet-1_pInt)
   endif
 enddo

640 do i = 1_pInt,elemSet
   if (mesh_mapElemSet(1,i) == 0_pInt) call IO_error(error_ID=904_pInt,ext_msg=mesh_nameElemSet(i))
 enddo

end subroutine mesh_abaqus_map_elementSets


!********************************************************************
! map solid section (Abaqus only)
!
! allocate globals: mesh_nameMaterial, mesh_mapMaterial
!********************************************************************
subroutine mesh_abaqus_map_materials(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 20_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt) :: i,c = 0_pInt
 logical :: inPart = .false.
 character(len=64) :: elemSetName,materialName
 
 allocate (mesh_nameMaterial(mesh_Nmaterials)) ; mesh_nameMaterial = ''
 allocate (mesh_mapMaterial(mesh_Nmaterials)) ;  mesh_mapMaterial = ''

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,myPos,1_pInt)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'section' ) then

     elemSetName = ''
     materialName = ''

     do i = 3_pInt,myPos(1_pInt)
       if (IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'elset') /= '') &
         elemSetName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'elset'))
       if (IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'material') /= '') &
         materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'material'))
     enddo

     if (elemSetName /= '' .and. materialName /= '') then
       c = c + 1_pInt
       mesh_nameMaterial(c) = materialName                                                          ! name of material used for this section
       mesh_mapMaterial(c)  = elemSetName                                                           ! mapped to respective element set
     endif       
   endif
 enddo

620 if (c==0_pInt) call IO_error(error_ID=905_pInt)
 do i=1_pInt,c
   if (mesh_nameMaterial(i)=='' .or. mesh_mapMaterial(i)=='') call IO_error(error_ID=905_pInt)
 enddo

 end subroutine mesh_abaqus_map_materials
 

!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of CP elements in mesh and stores them in 'mesh_NcpElems'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_cpElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_extractValue
                 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line
 integer(pInt) :: i,k
 logical :: materialFound = .false.
 character(len=64) ::materialName,elemSetName
 
 mesh_NcpElems = 0_pInt
 
610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_stringValue(line,myPos,1_pInt)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if (IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'material' .and. materialFound) then
         do i = 1_pInt,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) &                                            ! matched?
                 mesh_NcpElems = mesh_NcpElems + mesh_mapElemSet(1,k)                               ! add those elem count
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo
 
620 if (mesh_NcpElems == 0_pInt) call IO_error(error_ID=906_pInt)

end subroutine mesh_abaqus_count_cpElements


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPelem'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_elements(myUnit)

 use math, only: math_qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error
                 
 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) ::i,j,k,cpElem = 0_pInt
 logical :: materialFound = .false.
 character (len=64) materialName,elemSetName                                                        ! why limited to 64? ABAQUS?

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=660) line
   myPos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_stringValue(line,myPos,1_pInt)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'material' .and. materialFound) then
         do i = 1_pInt,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) then                                         ! matched?
                 do j = 1_pInt,mesh_mapElemSet(1,k)
                   cpElem = cpElem + 1_pInt
                   mesh_mapFEtoCPelem(1,cpElem) = mesh_mapElemSet(1_pInt+j,k)                       ! store FE id
                   mesh_mapFEtoCPelem(2,cpElem) = cpElem                                            ! store our id
                 enddo
               endif
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo

660 call math_qsort(mesh_mapFEtoCPelem,1_pInt,int(size(mesh_mapFEtoCPelem,2_pInt),pInt))                  ! should be mesh_NcpElems

 if (int(size(mesh_mapFEtoCPelem),pInt) < 2_pInt) call IO_error(error_ID=907_pInt)

end subroutine mesh_abaqus_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPnode'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_nodes(myUnit)

 use math, only: math_qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_intValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt) :: i,c,cpNode = 0_pInt
 logical :: inPart = .false.

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=650) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
   ) then
     c = IO_countDataLines(myUnit)
     do i = 1_pInt,c
       backspace(myUnit)
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=650) line
       myPos = IO_stringPos(line,maxNchunks)
       cpNode = cpNode + 1_pInt
       mesh_mapFEtoCPnode(1_pInt,cpNode) = IO_intValue(line,myPos,1_pInt)
       mesh_mapFEtoCPnode(2_pInt,cpNode) = cpNode
     enddo
   endif
 enddo

650 call math_qsort(mesh_mapFEtoCPnode,1_pInt,int(size(mesh_mapFEtoCPnode,2_pInt),pInt))

 if (int(size(mesh_mapFEtoCPnode),pInt) == 0_pInt) call IO_error(error_ID=908_pInt)

end subroutine mesh_abaqus_map_nodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!! Allocates global arrays 'mesh_node0' and 'mesh_node'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_build_nodes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_floatValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_countDataLines, &
                 IO_intValue
use numerics, only: numerics_unitlength

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,j,m,c
 logical :: inPart
 
 allocate ( mesh_node0 (3,mesh_Nnodes) ); mesh_node0 = 0.0_pReal
 allocate ( mesh_node  (3,mesh_Nnodes) ); mesh_node  = 0.0_pReal

610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do
   read (myUnit,610,END=670) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
   ) then
     c = IO_countDataLines(myUnit)                                                                  ! how many nodes are defined here?
     do i = 1_pInt,c
       backspace(myUnit)                                                                            ! rewind to first entry
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=670) line
       myPos = IO_stringPos(line,maxNchunks)
       m = mesh_FEasCP('node',IO_intValue(line,myPos,1_pInt))
       do j=1_pInt, 3_pInt
         mesh_node0(j,m) = numerics_unitlength * IO_floatValue(line,myPos,j+1_pInt)
       enddo  
     enddo
   endif
 enddo

670 if (int(size(mesh_node0,2_pInt),pInt) /= mesh_Nnodes) call IO_error(error_ID=909_pInt)
 mesh_node = mesh_node0

end subroutine mesh_abaqus_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets maximum count of nodes, IPs, IP neighbors, and subNodes among cpElements.
!! Allocates global arrays 'mesh_maxNnodes', 'mesh_maxNips', mesh_maxNipNeighbors', 
!! and mesh_maxNsubNodes
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_cpSizes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue ,&
                 IO_error, &
                 IO_countDataLines, &
                 IO_intValue

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,c,t,g
 logical :: inPart

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNsubNodes    = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'type'))           ! remember elem type
     g = FE_geomtype(t)
     mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(g))
     mesh_maxNips =         max(mesh_maxNips,FE_Nips(g))
     mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(g))
     mesh_maxNsubNodes =    max(mesh_maxNsubNodes,FE_NsubNodes(g))
   endif
 enddo
 
620 end subroutine mesh_abaqus_count_cpSizes


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, mat, tex, and node list per elemen.
!! Allocates global array 'mesh_element'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_build_elements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_skipChunks, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_extractValue, &
                 IO_floatValue, &
                 IO_error, &
                 IO_countDataLines

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 65_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos

 integer(pInt) :: i,j,k,c,e,t,homog,micro
 logical inPart,materialFound
 character (len=64) :: materialName,elemSetName
 character(len=300) :: line

 allocate (mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do
   read (myUnit,610,END=620) line
   myPos(1:1+2*2) = IO_stringPos(line,2_pInt)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'type'))          ! remember elem type
     c = IO_countDataLines(myUnit)
     do i = 1_pInt,c
       backspace(myUnit)
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=620) line
       myPos = IO_stringPos(line,maxNchunks)                                                       ! limit to 64 nodes max
       e = mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt))
       if (e /= 0_pInt) then                                                                       ! disregard non CP elems
         mesh_element(1,e) = IO_intValue(line,myPos,1_pInt)                                        ! FE id
         mesh_element(2,e) = t                                                                     ! elem type
         do j=1_pInt,FE_Nnodes(FE_geomtype(t))
           mesh_element(4_pInt+j,e) = IO_intValue(line,myPos,1_pInt+j)                             ! copy FE ids of nodes to position 5:
         enddo
         call IO_skipChunks(myUnit,FE_NoriginalNodes(t)-(myPos(1_pInt)-1_pInt))                    ! read on (even multiple lines) if FE_NoriginalNodes exceeds required node count
       endif
     enddo
   endif
 enddo

 
620 rewind(myUnit)                                                                                 ! just in case "*material" definitions apear before "*element"

 materialFound = .false.
 do 
   read (myUnit,610,END=630) line
   myPos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_StringValue(line,myPos,1_pInt)))
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_StringValue(line,myPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if ( IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'material' .and. &
            materialFound ) then
         read (myUnit,610,END=630) line                                                             ! read homogenization and microstructure
         myPos(1:1+2*2) = IO_stringPos(line,2_pInt)
         homog = nint(IO_floatValue(line,myPos,1_pInt),pInt)
         micro = nint(IO_floatValue(line,myPos,2_pInt),pInt)
         do i = 1_pInt,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) then                                         ! matched?
                 do j = 1_pInt,mesh_mapElemSet(1,k)
                   e = mesh_FEasCP('elem',mesh_mapElemSet(1+j,k))
                   mesh_element(3,e) = homog                                                        ! store homogenization
                   mesh_element(4,e) = micro                                                        ! store microstructure
                   mesh_maxValStateVar(1) = max(mesh_maxValStateVar(1),homog)
                   mesh_maxValStateVar(2) = max(mesh_maxValStateVar(2),micro)
                 enddo
               endif
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo

630 end subroutine mesh_abaqus_build_elements
#endif


!********************************************************************
! get any additional damask options from input file
!
! mesh_periodicSurface
!********************************************************************
subroutine mesh_get_damaskOptions(myUnit)

use IO, only: &
  IO_lc, &
  IO_stringValue, &
  IO_stringPos

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 5_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 integer(pInt) chunk, Nchunks
 character(len=300) :: line, damaskOption, v
#ifndef Spectral
 character(len=300) :: keyword 
#endif
 mesh_periodicSurface = .false.

610 FORMAT(A300)

#ifdef Marc 
 keyword = '$damask'
#endif
#ifdef Abaqus
 keyword = '**damask'
#endif

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   Nchunks = myPos(1)
#ifndef Spectral
   if (IO_lc(IO_stringValue(line,myPos,1_pInt)) == keyword .and. Nchunks > 1_pInt) then  ! found keyword for damask option and there is at least one more chunk to read
     damaskOption = IO_lc(IO_stringValue(line,myPos,2_pInt))
     select case(damaskOption)
       case('periodic')                                                        ! damask Option that allows to specify periodic fluxes
         do chunk = 3_pInt,Nchunks                                                  ! loop through chunks (skipping the keyword)
            v = IO_lc(IO_stringValue(line,myPos,chunk))                       ! chunk matches keyvalues x,y, or z?
            mesh_periodicSurface(1) = mesh_periodicSurface(1) .or. v == 'x'
            mesh_periodicSurface(2) = mesh_periodicSurface(2) .or. v == 'y'
            mesh_periodicSurface(3) = mesh_periodicSurface(3) .or. v == 'z'
         enddo
     endselect
   endif
#else
   damaskOption = IO_lc(IO_stringValue(line,myPos,1_pInt))
   select case(damaskOption)
     case('periodic')                                                        ! damask Option that allows to specify periodic fluxes
       do chunk = 2_pInt,Nchunks                                                  ! loop through chunks (skipping the keyword)
          v = IO_lc(IO_stringValue(line,myPos,chunk))                       ! chunk matches keyvalues x,y, or z?
          mesh_periodicSurface(1) = mesh_periodicSurface(1) .or. v == 'x'
          mesh_periodicSurface(2) = mesh_periodicSurface(2) .or. v == 'y'
          mesh_periodicSurface(3) = mesh_periodicSurface(3) .or. v == 'z'
        enddo
   endselect
#endif
 enddo

620 end subroutine mesh_get_damaskOptions


!***********************************************************
! calculation of IP interface areas
!
! allocate globals
! _ipArea, _ipAreaNormal
!***********************************************************
subroutine mesh_build_ipAreas
 
 use math, only: math_vectorproduct
 
 implicit none
 integer(pInt) :: e,f,t,i,j,n,Ntriangles                                     ! each interface is made up of this many triangles
 real(pReal), dimension (3,FE_maxNipFaceNodes) :: nPos                       ! coordinates of nodes on IP face
 real(pReal), dimension(3,FE_maxNipFaceNodes-2_pInt,FE_maxNipFaceNodes) :: normal
 real(pReal), dimension(FE_maxNipFaceNodes-2_pInt,FE_maxNipFaceNodes) :: area

 allocate(mesh_ipArea(mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ;              mesh_ipArea       = 0.0_pReal
 allocate(mesh_ipAreaNormal(3_pInt,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ; mesh_ipAreaNormal = 0.0_pReal
 do e = 1_pInt,mesh_NcpElems                                                 ! loop over cpElems
   t = FE_geomtype(mesh_element(2,e))                                        ! get elemGeomType
   select case (FE_dimension(t))                                             ! treat (1D), 2D, and 3D element types differently

     case (2)
       do i = 1_pInt,FE_Nips(t)                                              ! loop over IPs of elem
         do f = 1_pInt,FE_NipNeighbors(t)                                    ! loop over interfaces of IP
           forall (n = 1_pInt:2_pInt) nPos(1:3,n) = mesh_subNodeCoord(1:3,FE_subNodeOnIPFace(n,f,i,t),e)
           mesh_ipAreaNormal(1,f,i,e) =   nPos(2,2) - nPos(2,1)                                                      ! x_normal =  y_connectingVector
           mesh_ipAreaNormal(2,f,i,e) = -(nPos(1,2) - nPos(1,1))                                                     ! y_normal = -x_connectingVector
           mesh_ipArea(f,i,e) = sqrt(sum(mesh_ipAreaNormal(1:3,f,i,e)*mesh_ipAreaNormal(1:3,f,i,e)))                 ! and area
           mesh_ipAreaNormal(1:3,f,i,e) = mesh_ipAreaNormal(1:3,f,i,e) / mesh_ipArea(f,i,e)                          ! have unit normal
         enddo
       enddo
           
     case (3)
       Ntriangles = FE_NipFaceNodes(t) - 2_pInt
       do i = 1_pInt,FE_Nips(t)                                             ! loop over IPs of elem
         do f = 1_pInt,FE_NipNeighbors(t)                                   ! loop over interfaces of IP
           nPos = 0.0_pReal
           normal = 0.0_pReal
           area = 0.0_pReal
           forall (n = 1_pInt:FE_NipFaceNodes(t)) nPos(1:3,n) = mesh_subNodeCoord(1:3,FE_subNodeOnIPFace(n,f,i,t),e)
           forall (n = 1_pInt:FE_NipFaceNodes(t), j = 1_pInt:Ntriangles)    ! start at each interface node and build valid triangles to cover interface
             normal(1:3,j,n) = math_vectorproduct(nPos(1:3,1_pInt+mod(n+j-1_pInt,FE_NipFaceNodes(t))) - nPos(1:3,n), &    ! calc their normal vectors
                                                  nPos(1:3,1_pInt+mod(n+j-0_pInt,FE_NipFaceNodes(t))) - nPos(1:3,n))
             area(j,n) = sqrt(sum(normal(1:3,j,n)*normal(1:3,j,n)))                                                    ! and area
           end forall
           forall (n = 1_pInt:FE_NipFaceNodes(t), j = 1_pInt:Ntriangles, area(j,n) > 0.0_pReal) &
             normal(1:3,j,n) = normal(1:3,j,n) / area(j,n)                                     ! have each normal of unit length
           mesh_ipArea(f,i,e) = sum(area) / (FE_NipFaceNodes(t)*2.0_pReal)                        ! vector product gives area of parallelograms instead of triangles
           mesh_ipAreaNormal(1:3,f,i,e) = sum(sum(normal,3),2)/ &
                                          max(1.0_pReal,real(count(area > 0.0_pReal),pReal))   ! average of all valid normals (not necessarily unit length---beware!!/fix?)
         enddo
       enddo

    end select
  enddo
 
 end subroutine mesh_build_ipAreas
 
 
!***********************************************************
! assignment of twin nodes for each cp node
!
! allocate globals
! _nodeTwins
!***********************************************************
subroutine mesh_build_nodeTwins

implicit none
integer(pInt) dir, &      ! direction of periodicity
              node, &
              minimumNode, &
              maximumNode, &
              n1, &
              n2
integer(pInt), dimension(mesh_Nnodes+1) :: minimumNodes, maximumNodes ! list of surface nodes (minimum and maximum coordinate value) with first entry giving the number of nodes
real(pReal)   minCoord, maxCoord, &     ! extreme positions in one dimension
              tolerance                 ! tolerance below which positions are assumed identical
real(pReal), dimension(3) ::  distance  ! distance between two nodes in all three coordinates
logical, dimension(mesh_Nnodes) :: unpaired

allocate(mesh_nodeTwins(3,mesh_Nnodes))
mesh_nodeTwins = 0_pInt

tolerance = 0.001_pReal * minval(mesh_ipVolume) ** 0.333_pReal

do dir = 1_pInt,3_pInt                                    ! check periodicity in directions of x,y,z
  if (mesh_periodicSurface(dir)) then                     ! only if periodicity is requested

    
    !*** find out which nodes sit on the surface 
    !*** and have a minimum or maximum position in this dimension
    
    minimumNodes = 0_pInt
    maximumNodes = 0_pInt
    minCoord = minval(mesh_node0(dir,:))
    maxCoord = maxval(mesh_node0(dir,:))
    do node = 1_pInt,mesh_Nnodes                     ! loop through all nodes and find surface nodes
      if (abs(mesh_node0(dir,node) - minCoord) <= tolerance) then
        minimumNodes(1) = minimumNodes(1) + 1_pInt
        minimumNodes(minimumNodes(1)+1_pInt) = node
      elseif (abs(mesh_node0(dir,node) - maxCoord) <= tolerance) then
        maximumNodes(1) = maximumNodes(1) + 1_pInt
        maximumNodes(maximumNodes(1)+1_pInt) = node
      endif
    enddo
    
    
    !*** find the corresponding node on the other side with the same position in this dimension
    
    unpaired = .true.
    do n1 = 1_pInt,minimumNodes(1)
      minimumNode = minimumNodes(n1+1_pInt)
      if (unpaired(minimumNode)) then
        do n2 = 1_pInt,maximumNodes(1)
          maximumNode = maximumNodes(n2+1_pInt)
          distance = abs(mesh_node0(:,minimumNode) - mesh_node0(:,maximumNode))
          if (sum(distance) - distance(dir) <= tolerance) then        ! minimum possible distance (within tolerance)
            mesh_nodeTwins(dir,minimumNode) = maximumNode
            mesh_nodeTwins(dir,maximumNode) = minimumNode
            unpaired(maximumNode) = .false.                           ! remember this node, we don't have to look for his partner again
            exit
          endif
        enddo
      endif
    enddo

  endif
enddo

end subroutine mesh_build_nodeTwins



 






!********************************************************************
! get maximum count of shared elements among cpElements and
! build list of elements shared by each node in mesh
!
! _maxNsharedElems
! _sharedElem
!********************************************************************
subroutine mesh_build_sharedElems

implicit none
integer(pint)   e, &                                                                                ! element index
                t, &                                                                                ! element type
                node, &                                                                             ! CP node index
                j, &                                                                                ! node index per element 
                myDim, &                                                                            ! dimension index 
                nodeTwin                                                                            ! node twin in the specified dimension
integer(pInt), dimension (mesh_Nnodes) :: node_count
integer(pInt), dimension (:), allocatable :: node_seen

allocate(node_seen(maxval(FE_Nnodes)))


node_count = 0_pInt

do e = 1_pInt,mesh_NcpElems
  t = FE_geomtype(mesh_element(2,e))                                                                ! get elemGeomType
  node_seen = 0_pInt                                                                                ! reset node duplicates
  do j = 1_pInt,FE_Nnodes(t)                                                                             ! check each node of element
    node = mesh_FEasCP('node',mesh_element(4+j,e))                                                  ! translate to internal (consecutive) numbering
    if (all(node_seen /= node)) then
      node_count(node) = node_count(node) + 1_pInt                                                  ! if FE node not yet encountered -> count it
      do myDim = 1_pInt,3_pInt                                                                                  ! check in each dimension...
        nodeTwin = mesh_nodeTwins(myDim,node)
        if (nodeTwin > 0_pInt) &                                                                    ! if I am a twin of some node...
          node_count(nodeTwin) = node_count(nodeTwin) + 1_pInt                                      ! -> count me again for the twin node
      enddo
    endif
    node_seen(j) = node                                                                             ! remember this node to be counted already
  enddo
enddo

mesh_maxNsharedElems = int(maxval(node_count),pInt)                                                 ! most shared node

allocate(mesh_sharedElem(1+mesh_maxNsharedElems,mesh_Nnodes))
mesh_sharedElem = 0_pInt

do e = 1_pInt,mesh_NcpElems
  t = FE_geomtype(mesh_element(2,e))                                                                ! get elemGeomType
  node_seen = 0_pInt
  do j = 1_pInt,FE_Nnodes(t)
    node = mesh_FEasCP('node',mesh_element(4_pInt+j,e))
    if (all(node_seen /= node)) then
      mesh_sharedElem(1,node) = mesh_sharedElem(1,node) + 1_pInt                                    ! count for each node the connected elements
      mesh_sharedElem(mesh_sharedElem(1,node)+1_pInt,node) = e                                      ! store the respective element id
      do myDim = 1_pInt,3_pInt                                                                      ! check in each dimension...
        nodeTwin = mesh_nodeTwins(myDim,node)
        if (nodeTwin > 0_pInt) then                                                                 ! if i am a twin of some node...
          mesh_sharedElem(1,nodeTwin) = mesh_sharedElem(1,nodeTwin) + 1_pInt                        ! ...count me again for the twin
          mesh_sharedElem(mesh_sharedElem(1,nodeTwin)+1,nodeTwin) = e                               ! store the respective element id
        endif
      enddo
    endif
    node_seen(j) = node
  enddo
enddo

deallocate(node_seen)

end subroutine mesh_build_sharedElems


!***********************************************************
! build up of IP neighborhood
!
! allocate globals
! _ipNeighborhood
!***********************************************************
subroutine mesh_build_ipNeighborhood

use math, only: math_mul3x3

implicit none
integer(pInt)                   myElem, &                                                           ! my CP element index
                                myIP, &
                                myType, &                                                           ! my element type
                                myFace, &
                                neighbor, &                                                         ! neighor index
                                neighboringIPkey, &                                                 ! positive integer indicating the neighboring IP (for intra-element) and negative integer indicating the face towards neighbor (for neighboring element) 
                                candidateIP, &
                                neighboringType, &                                                  ! element type of neighbor
                                NlinkedNodes, &                                                     ! number of linked nodes
                                twin_of_linkedNode, &                                               ! node twin of a specific linkedNode
                                NmatchingNodes, &                                                   ! number of matching nodes
                                dir, &                                                              ! direction of periodicity
                                matchingElem, &                                                     ! CP elem number of matching element
                                matchingFace, &                                                     ! face ID of matching element
                                a, anchor, &
                                neighboringIP, &  
                                neighboringElem, &
                                pointingToMe
integer(pInt), dimension(FE_maxmaxNnodesAtIP) :: &
                                linkedNodes = 0_pInt, &
                                matchingNodes
logical checkTwins

allocate(mesh_ipNeighborhood(3,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems))
mesh_ipNeighborhood = 0_pInt


do myElem = 1_pInt,mesh_NcpElems                                                                    ! loop over cpElems
  myType = FE_geomtype(mesh_element(2,myElem))                                                      ! get elemGeomType
  do myIP = 1_pInt,FE_Nips(myType)                                                                  ! loop over IPs of elem

    do neighbor = 1_pInt,FE_NipNeighbors(myType)                                                    ! loop over neighbors of IP
      neighboringIPkey = FE_ipNeighbor(neighbor,myIP,myType)

      !*** if the key is positive, the neighbor is inside the element
      !*** that means, we have already found our neighboring IP
      
      if (neighboringIPkey > 0_pInt) then
        mesh_ipNeighborhood(1,neighbor,myIP,myElem) = myElem
        mesh_ipNeighborhood(2,neighbor,myIP,myElem) = neighboringIPkey


      !*** if the key is negative, the neighbor resides in a neighboring element
      !*** that means, we have to look through the face indicated by the key and see which element is behind that face
      
      elseif (neighboringIPkey < 0_pInt) then                                            ! neighboring element's IP
        myFace = -neighboringIPkey
        call mesh_faceMatch(myElem, myFace, matchingElem, matchingFace)                  ! get face and CP elem id of face match
        if (matchingElem > 0_pInt) then                                                  ! found match?
          neighboringType = FE_geomtype(mesh_element(2,matchingElem))

          !*** trivial solution if neighbor has only one IP
          
          if (FE_Nips(neighboringType) == 1_pInt) then            
            mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
            mesh_ipNeighborhood(2,neighbor,myIP,myElem) = 1_pInt
            cycle
          endif

          !*** find those nodes which build the link to the neighbor
          
          NlinkedNodes = 0_pInt
          linkedNodes = 0_pInt
          do a = 1_pInt,FE_maxNnodesAtIP(myType)                                         ! figure my anchor nodes on connecting face
            anchor = FE_nodesAtIP(a,myIP,myType)
            if (anchor /= 0_pInt) then                                                   ! valid anchor node
              if (any(FE_nodeOnFace(:,myFace,myType) == anchor)) then                    ! ip anchor sits on face?
                NlinkedNodes = NlinkedNodes + 1_pInt
                linkedNodes(NlinkedNodes) = &
                   mesh_FEasCP('node',mesh_element(4_pInt+anchor,myElem))                ! CP id of anchor node
              else                                                                       ! something went wrong with the linkage, since not all anchors sit on my face
                NlinkedNodes = 0_pInt
                linkedNodes = 0_pInt
                exit
              endif
            endif
          enddo

          !*** loop through the ips of my neighbor
          !*** and try to find an ip with matching nodes
          !*** also try to match with node twins

checkCandidateIP: do candidateIP = 1_pInt,FE_Nips(neighboringType)
            NmatchingNodes = 0_pInt
            matchingNodes = 0_pInt
            do a = 1_pInt,FE_maxNnodesAtIP(neighboringType)                              ! check each anchor node of that ip
              anchor = FE_nodesAtIP(a,candidateIP,neighboringType)
              if (anchor /= 0_pInt) then                                            ! valid anchor node
                if (any(FE_nodeOnFace(:,matchingFace,neighboringType) == anchor)) then  ! sits on matching face?
                  NmatchingNodes = NmatchingNodes + 1_pInt
                  matchingNodes(NmatchingNodes) = &
                     mesh_FEasCP('node',mesh_element(4+anchor,matchingElem))        ! CP id of neighbor's anchor node
                else                                                                ! no matching, because not all nodes sit on the matching face
                  NmatchingNodes = 0_pInt
                  matchingNodes = 0_pInt
                  exit
                endif
              endif
            enddo

            if (NmatchingNodes /= NlinkedNodes) &                                   ! this ip has wrong count of anchors on face
              cycle checkCandidateIP
            
            !*** check "normal" nodes whether they match or not
            
            checkTwins = .false.
            do a = 1_pInt,NlinkedNodes
              if (all(matchingNodes /= linkedNodes(a))) then                        ! this linkedNode does not match any matchingNode
                checkTwins = .true.
                exit                                                                ! no need to search further
              endif
            enddo
            
            !*** if no match found, then also check node twins
            
            if(checkTwins) then
              dir = int(maxloc(abs(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem)),1),pInt)      ! check for twins only in direction of the surface normal
              do a = 1_pInt,NlinkedNodes
                twin_of_linkedNode = mesh_nodeTwins(dir,linkedNodes(a))
                if (twin_of_linkedNode == 0_pInt .or. &                             ! twin of linkedNode does not exist...
                    all(matchingNodes /= twin_of_linkedNode)) then                  ! ... or it does not match any matchingNode
                  cycle checkCandidateIP                                            ! ... then check next candidateIP
                endif
              enddo
            endif

            !*** we found a match !!!

            mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
            mesh_ipNeighborhood(2,neighbor,myIP,myElem) = candidateIP
            exit checkCandidateIP            
          enddo checkCandidateIP
        endif                                                                       ! end of valid external matching
      endif                                                                         ! end of internal/external matching
    enddo
  enddo
enddo
do myElem = 1_pInt,mesh_NcpElems                                                                                 ! loop over cpElems
  myType = FE_geomtype(mesh_element(2,myElem))                                                                   ! get elemGeomType
  do myIP = 1_pInt,FE_Nips(myType)                                                                               ! loop over IPs of elem
    do neighbor = 1_pInt,FE_NipNeighbors(myType)                                                                 ! loop over neighbors of IP
      neighboringElem = mesh_ipNeighborhood(1,neighbor,myIP,myElem)
      neighboringIP   = mesh_ipNeighborhood(2,neighbor,myIP,myElem)
      if (neighboringElem > 0_pInt .and. neighboringIP > 0_pInt) then                                            ! if neighbor exists ...
        neighboringType = FE_geomtype(mesh_element(2,neighboringElem))
        do pointingToMe = 1_pInt,FE_NipNeighbors(neighboringType)                                                ! find neighboring index that points from my neighbor to myself
          if (    myElem == mesh_ipNeighborhood(1,pointingToMe,neighboringIP,neighboringElem) &
              .and. myIP == mesh_ipNeighborhood(2,pointingToMe,neighboringIP,neighboringElem)) then              ! possible candidate
            if (math_mul3x3(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem),&
                            mesh_ipAreaNormal(1:3,pointingToMe,neighboringIP,neighboringElem)) < 0.0_pReal) then ! area normals have opposite orientation (we have to check that because of special case for single element with two ips and periodicity. In this case the neighbor is identical in two different directions.)
              mesh_ipNeighborhood(3,neighbor,myIP,myElem) = pointingToMe                                         ! found match
              exit                                                                                               ! so no need to search further
            endif
          endif
        enddo
      endif
    enddo
  enddo
enddo

end subroutine mesh_build_ipNeighborhood


!***********************************************************
! write statistics regarding input file parsing
! to the output file
! 
!***********************************************************
subroutine mesh_tell_statistics

 use math,  only: math_range
 use IO,    only: IO_error
 use debug, only: debug_level, &
                  debug_mesh, &
                  debug_levelBasic, &
                  debug_levelExtensive, &
                  debug_levelSelective, &
                  debug_e, &
                  debug_i

 implicit none
 integer(pInt), dimension (:,:), allocatable :: mesh_HomogMicro
 character(len=64) :: myFmt
 integer(pInt) :: i,e,n,f,t, myDebug
 
 myDebug = debug_level(debug_mesh)

 if (mesh_maxValStateVar(1) < 1_pInt) call IO_error(error_ID=170_pInt) ! no homogenization specified
 if (mesh_maxValStateVar(2) < 1_pInt) call IO_error(error_ID=180_pInt) ! no microstructure specified
 
 allocate (mesh_HomogMicro(mesh_maxValStateVar(1),mesh_maxValStateVar(2))); mesh_HomogMicro = 0_pInt
do e = 1_pInt,mesh_NcpElems
  if (mesh_element(3,e) < 1_pInt) call IO_error(error_ID=170_pInt,e=e) ! no homogenization specified
  if (mesh_element(4,e) < 1_pInt) call IO_error(error_ID=180_pInt,e=e) ! no microstructure specified
  mesh_HomogMicro(mesh_element(3,e),mesh_element(4,e)) = &
  mesh_HomogMicro(mesh_element(3,e),mesh_element(4,e)) + 1_pInt ! count combinations of homogenization and microstructure
enddo
!$OMP CRITICAL (write2out)
  if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
    write(6,'(/,a,/)') ' Input Parser: STATISTICS'
    write(6,*) mesh_Nelems,           ' : total number of elements in mesh'
    write(6,*) mesh_NcpElems,         ' : total number of CP elements in mesh'
    write(6,*) mesh_Nnodes,           ' : total number of nodes in mesh'
    write(6,*) mesh_maxNnodes,        ' : max number of nodes in any CP element'
    write(6,*) mesh_maxNips,          ' : max number of IPs in any CP element'
    write(6,*) mesh_maxNipNeighbors,  ' : max number of IP neighbors in any CP element'
    write(6,*) mesh_maxNsubNodes,     ' : max number of (additional) subnodes in any CP element'
    write(6,*) mesh_maxNsharedElems,  ' : max number of CP elements sharing a node'
    write(6,'(/,a,/)') ' Input Parser: HOMOGENIZATION/MICROSTRUCTURE'
    write(6,*) mesh_maxValStateVar(1), ' : maximum homogenization index'
    write(6,*) mesh_maxValStateVar(2), ' : maximum microstructure index'
    write(6,*)
    write (myFmt,'(a,i32.32,a)') '(9x,a2,1x,',mesh_maxValStateVar(2),'(i8))'
    write(6,myFmt) '+-',math_range(mesh_maxValStateVar(2))
    write (myFmt,'(a,i32.32,a)') '(i8,1x,a2,1x,',mesh_maxValStateVar(2),'(i8))'
    do i=1_pInt,mesh_maxValStateVar(1)      ! loop over all (possibly assigned) homogenizations
      write(6,myFmt) i,'| ',mesh_HomogMicro(i,:) ! loop over all (possibly assigned) microstructures
    enddo
    write(6,'(/,a,/)') ' Input Parser: ADDITIONAL MPIE OPTIONS'
    write(6,*) 'periodic surface : ', mesh_periodicSurface
    write(6,*)
    flush(6)
  endif

  if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
    write(6,*)
    write(6,*) 'Input Parser: SUBNODE COORDINATES'
    write(6,*)
    write(6,'(a8,1x,a5,1x,2(a15,1x),a20,3(1x,a12))')&
                              'elem','IP','IP neighbor','IPFaceNodes','subNodeOnIPFace','x','y','z'
    do e = 1_pInt,mesh_NcpElems                                                          ! loop over cpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = FE_geomtype(mesh_element(2,e))                                                 ! get elemGeomType
      do i = 1_pInt,FE_Nips(t)                                                           ! loop over IPs of elem
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        do f = 1_pInt,FE_NipNeighbors(t)                                                 ! loop over interfaces of IP
          do n = 1_pInt,FE_NipFaceNodes(t)                                                  ! loop over nodes on interface
            write(6,'(i8,1x,i5,2(1x,i15),1x,i20,3(1x,f12.8))') e,i,f,n,FE_subNodeOnIPFace(n,f,i,t),&
                                               mesh_subNodeCoord(1,FE_subNodeOnIPFace(n,f,i,t),e),&
                                               mesh_subNodeCoord(2,FE_subNodeOnIPFace(n,f,i,t),e),&
                                               mesh_subNodeCoord(3,FE_subNodeOnIPFace(n,f,i,t),e)
          enddo
        enddo
      enddo
    enddo
    write(6,*)
    write(6,*) 'Input Parser: IP COORDINATES'
    write(6,'(a8,1x,a5,3(1x,a12))') 'elem','IP','x','y','z'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write(6,'(i8,1x,i5,3(1x,f12.8))') e, i, mesh_ipCoordinates(:,i,e)
      enddo
    enddo 
    write(6,*)
    write(6,*) 'Input Parser: ELEMENT VOLUME'
    write(6,*)
    write(6,'(a13,1x,e15.8)') 'total volume', sum(mesh_ipVolume)
    write(6,*)
    write(6,'(a8,1x,a5,1x,a15,1x,a5,1x,a15,1x,a16)') 'elem','IP','volume','face','area','-- normal --'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = FE_geomtype(mesh_element(2,e))         ! get elemGeomType
      do i = 1_pInt,FE_Nips(t)
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write(6,'(i8,1x,i5,1x,e15.8)') e,i,mesh_IPvolume(i,e)
        do f = 1_pInt,FE_NipNeighbors(t)
          write(6,'(i33,1x,e15.8,1x,3(f6.3,1x))') f,mesh_ipArea(f,i,e),mesh_ipAreaNormal(:,f,i,e)
        enddo
      enddo
    enddo
    write(6,*)
    write(6,*) 'Input Parser: NODE TWINS'
    write(6,*)
    write(6,'(a6,3(3x,a6))') '  node','twin_x','twin_y','twin_z'
    do n = 1_pInt,mesh_Nnodes                    ! loop over cpNodes
      if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. .not. any(mesh_element(5:,debug_e) == n)) cycle
      write(6,'(i6,3(3x,i6))') n, mesh_nodeTwins(1:3,n)
    enddo
    write(6,*)
    write(6,*) 'Input Parser: IP NEIGHBORHOOD'
    write(6,*)
    write(6,'(a8,1x,a10,1x,a10,1x,a3,1x,a13,1x,a13)') 'elem','IP','neighbor','','elemNeighbor','ipNeighbor'
    do e = 1_pInt,mesh_NcpElems                                                          ! loop over cpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = FE_geomtype(mesh_element(2,e))                                                 ! get elemGeomType
      do i = 1_pInt,FE_Nips(t)                                                           ! loop over IPs of elem
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        do n = 1_pInt,FE_NipNeighbors(t)                                                 ! loop over neighbors of IP
          write(6,'(i8,1x,i10,1x,i10,1x,a3,1x,i13,1x,i13)') e,i,n,'-->',mesh_ipNeighborhood(1,n,i,e),mesh_ipNeighborhood(2,n,i,e)
        enddo
      enddo
    enddo
  endif
!$OMP END CRITICAL (write2out)

 deallocate(mesh_HomogMicro)
 
end subroutine mesh_tell_statistics


!***********************************************************
! mapping of FE element types to internal representation
!***********************************************************
integer(pInt) function FE_mapElemtype(what)
 
 use IO, only: IO_lc, IO_error

 implicit none
 character(len=*), intent(in) :: what
  
 select case (IO_lc(what))
    case (   '6')
      FE_mapElemtype = 1_pInt            ! Two-dimensional Plane Strain Triangle
    case ( '155', &
           '125', &
           '128')
      FE_mapElemtype = 2_pInt            ! Two-dimensional Plane Strain triangle (155: cubic shape function, 125/128: second order isoparametric)
    case ( '11', &
           'cpe4')
      FE_mapElemtype = 3_pInt            ! Arbitrary Quadrilateral Plane-strain
    case ( '27', &
           'cpe8')
      FE_mapElemtype = 4_pInt            ! Plane Strain, Eight-node Distorted Quadrilateral
    case ( '54')
      FE_mapElemtype = 5_pInt            ! Plane Strain, Eight-node Distorted Quadrilateral with reduced integration
    case ('134', &
          'c3d4')
      FE_mapElemtype = 6_pInt            ! Three-dimensional Four-node Tetrahedron
    case ('157')
      FE_mapElemtype = 7_pInt            ! Three-dimensional, Low-order, Tetrahedron, Herrmann Formulations
    case ('127')
      FE_mapElemtype = 8_pInt            ! Three-dimensional Ten-node Tetrahedron
    case ('136', &
          'c3d6')
      FE_mapElemtype = 9_pInt            ! Three-dimensional Arbitrarily Distorted Pentahedral
    case ( '117', &
           '123', &
           'c3d8r')
      FE_mapElemtype = 10_pInt           ! Three-dimensional Arbitrarily Distorted linear hexahedral with reduced integration
    case (  '7', &
            'c3d8')
      FE_mapElemtype = 11_pInt           ! Three-dimensional Arbitrarily Distorted Brick
    case ( '57', &
           'c3d20r')
      FE_mapElemtype = 12_pInt           ! Three-dimensional Arbitrarily Distorted quad hexahedral with reduced integration
    case ( '21', &
           'c3d20')
      FE_mapElemtype = 13_pInt           ! Three-dimensional Arbitrarily Distorted quadratic hexahedral
    case default 
      call IO_error(error_ID=190_pInt,ext_msg=IO_lc(what))
 end select

end function FE_mapElemtype


!***********************************************************
! find face-matching element of same type
!***********************************************************
subroutine mesh_faceMatch(elem, face ,matchingElem, matchingFace)

implicit none
!*** output variables
integer(pInt), intent(out) ::     matchingElem, &                                                   ! matching CP element ID
                                  matchingFace                                                      ! matching FE face ID 

!*** input variables
integer(pInt), intent(in) ::      face, &                                                           ! FE face ID
                                  elem                                                              ! FE elem ID

!*** local variables
integer(pInt), dimension(FE_NfaceNodes(face,FE_geomtype(mesh_element(2,elem)))) :: &
                                  myFaceNodes                                                       ! global node ids on my face
integer(pInt)        ::           myType, &
                                  candidateType, &
                                  candidateElem, &
                                  candidateFace, &
                                  candidateFaceNode, &
                                  minNsharedElems, &
                                  NsharedElems, &
                                  lonelyNode = 0_pInt, &
                                  i, &
                                  n, &
                                  dir                                                               ! periodicity direction
integer(pInt), dimension(:), allocatable :: element_seen
logical checkTwins

matchingElem = 0_pInt
matchingFace = 0_pInt
minNsharedElems = mesh_maxNsharedElems + 1_pInt                                                     ! init to worst case
myType = FE_geomtype(mesh_element(2_pInt,elem))                                                     ! figure elemGeomType

do n = 1_pInt,FE_NfaceNodes(face,myType)                                                            ! loop over nodes on face
  myFaceNodes(n) = mesh_FEasCP('node',mesh_element(4_pInt+FE_nodeOnFace(n,face,myType),elem))       ! CP id of face node
  NsharedElems = mesh_sharedElem(1_pInt,myFaceNodes(n))                                             ! figure # shared elements for this node
  if (NsharedElems < minNsharedElems) then
    minNsharedElems = NsharedElems                                                                  ! remember min # shared elems
    lonelyNode = n                                                                                  ! remember most lonely node
  endif
enddo

allocate(element_seen(minNsharedElems))
element_seen = 0_pInt

checkCandidate: do i = 1_pInt,minNsharedElems                                                       ! iterate over lonelyNode's shared elements
  candidateElem = mesh_sharedElem(1_pInt+i,myFaceNodes(lonelyNode))                                 ! present candidate elem
  if (all(element_seen /= candidateElem)) then                                                      ! element seen for the first time?
    element_seen(i) = candidateElem
    candidateType = FE_geomtype(mesh_element(2_pInt,candidateElem))                                 ! figure elemGeomType of candidate
checkCandidateFace: do candidateFace = 1_pInt,FE_maxNipNeighbors                                    ! check each face of candidate
      if (FE_NfaceNodes(candidateFace,candidateType) /= FE_NfaceNodes(face,myType) &                ! incompatible face
          .or. (candidateElem == elem .and. candidateFace == face)) then                            ! this is my face
        cycle checkCandidateFace
      endif
      checkTwins = .false.
      do n = 1_pInt,FE_NfaceNodes(candidateFace,candidateType)                                      ! loop through nodes on face
        candidateFaceNode = mesh_FEasCP('node', mesh_element(4_pInt+FE_nodeOnFace(n,candidateFace,candidateType),candidateElem))
        if (all(myFaceNodes /= candidateFaceNode)) then             ! candidate node does not match any of my face nodes
          checkTwins = .true.                                       ! perhaps the twin nodes do match
          exit
        endif
      enddo
      if(checkTwins) then
checkCandidateFaceTwins: do dir = 1_pInt,3_pInt
          do n = 1_pInt,FE_NfaceNodes(candidateFace,candidateType)       ! loop through nodes on face
            candidateFaceNode = mesh_FEasCP('node', mesh_element(4+FE_nodeOnFace(n,candidateFace,candidateType),candidateElem))
            if (all(myFaceNodes /= mesh_nodeTwins(dir,candidateFaceNode))) then ! node twin does not match either
              if (dir == 3_pInt) then
                cycle checkCandidateFace
              else
                cycle checkCandidateFaceTwins                       ! try twins in next dimension
              endif
            endif
          enddo
          exit checkCandidateFaceTwins
        enddo checkCandidateFaceTwins
      endif
      matchingFace = candidateFace
      matchingElem = candidateElem
      exit checkCandidate                                           ! found my matching candidate
    enddo checkCandidateFace
  endif
enddo checkCandidate

deallocate(element_seen)

end subroutine mesh_faceMatch


!********************************************************************
! get properties of different types of finite elements
!
! assign globals:
! FE_nodesAtIP, FE_ipNeighbor, FE_subNodeParent, FE_subNodeOnIPFace
!********************************************************************
subroutine mesh_build_FEdata

 implicit none
 integer(pInt) :: me
 allocate(FE_nodesAtIP(FE_maxmaxNnodesAtIP,FE_maxNips,FE_Ngeomtypes)) ; FE_nodesAtIP = 0_pInt
 allocate(FE_ipNeighbor(FE_maxNipNeighbors,FE_maxNips,FE_Ngeomtypes)) ; FE_ipNeighbor = 0_pInt
 allocate(FE_subNodeParent(FE_maxNips,FE_maxNsubNodes,FE_Ngeomtypes)) ; FE_subNodeParent = 0_pInt
 allocate(FE_subNodeOnIPFace(FE_maxNipFaceNodes,FE_maxNipNeighbors,FE_maxNips,FE_Ngeomtypes)) ; FE_subNodeOnIPFace = 0_pInt
 
 ! fill FE_nodesAtIP with data
 me = 0_pInt

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    1,2,3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
    1,  &
    2,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  27 (2D 8node 9ip)
    reshape(int([&
    1,0,  &
    1,2,  &
    2,0,  &
    1,4,  &
    0,0,  &
    2,3,  &
    4,0,  &
    3,4,  &
    3,0   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
    1,2,3,4   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4,  &
    5,  &
    6   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
    1,2,3,4,5,6,7,8   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3,  &
    5,  &
    6,  &
    8,  &
    7   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  21 (3D 20node 27ip)
    reshape(int([&
    1,0, 0,0,  &
    1,2, 0,0,  &
    2,0, 0,0,  &
    1,4, 0,0,  &
    1,3, 2,4,  &
    2,3, 0,0,  &
    4,0, 0,0,  &
    3,4, 0,0,  &
    3,0, 0,0,  &
    1,5, 0,0,  &
    1,6, 2,5,  &
    2,6, 0,0,  &
    1,8, 4,5,  &
    0,0, 0,0,  &
    2,7, 3,6,  &
    4,8, 0,0,  &
    3,8, 4,7,  &
    3,7, 0,0,  &
    5,0, 0,0,  &
    5,6, 0,0,  &
    6,0, 0,0,  &
    5,8, 0,0,  &
    5,7, 6,8,  &
    6,7, 0,0,  &
    8,0, 0,0,  &
    7,8, 0,0,  &
    7,0, 0,0   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 
 ! *** FE_ipNeighbor ***
 ! is a list of the neighborhood of each IP.
 ! It is sorted in (local) +x,-x, +y,-y, +z,-z direction.
 ! Positive integers denote an intra-FE IP identifier.
 ! Negative integers denote the interface behind which the neighboring (extra-FE) IP will be located.
 me = 0_pInt

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    -2,-3,-1   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])
 
 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
     2,-3, 3,-1,  &
    -2, 1, 3,-1,  &
     2,-3,-2, 1   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])
 
 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
     2,-4, 3,-1,  &
    -2, 1, 4,-1,  &
     4,-4,-3, 1,  &
    -2, 3,-3, 2   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element  27 (2D 8node 9ip)
    reshape(int([&
     2,-4, 4,-1,  &
     3, 1, 5,-1,  &
    -2, 2, 6,-1,  &
     5,-4, 7, 1,  &
     6, 4, 8, 2,  &
    -2, 5, 9, 3,  &
     8,-4,-3, 4,  &
     9, 7,-3, 5,  &
    -2, 8,-3, 6   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
    -1,-2,-3,-4   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -2, 1, 3,-2, 4,-1,  &
     2,-4,-3, 1, 4,-1,  &
     2,-4, 3,-2,-3, 1   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -3, 1, 3,-2, 5,-1,  &
     2,-4,-3, 1, 6,-1,  &
     5,-4, 6,-2,-5, 1,  &
    -3, 4, 6,-2,-5, 2,  &
     5,-4,-3, 4,-5, 3   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
    -3,-5,-4,-2,-6,-1   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
     2,-5, 3,-2, 5,-1,  &
    -3, 1, 4,-2, 6,-1,  &
     4,-5,-4, 1, 7,-1,  &
    -3, 3,-4, 2, 8,-1,  &
     6,-5, 7,-2,-6, 1,  &
    -3, 5, 8,-2,-6, 2,  &
     8,-5,-4, 5,-6, 3,  &
    -3, 7,-4, 6,-6, 4   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element  21 (3D 20node 27ip)
    reshape(int([&
     2,-5, 4,-2,10,-1,  &
     3, 1, 5,-2,11,-1,  &
    -3, 2, 6,-2,12,-1,  &
     5,-5, 7, 1,13,-1,  &
     6, 4, 8, 2,14,-1,  &
    -3, 5, 9, 3,15,-1,  &
     8,-5,-4, 4,16,-1,  &
     9, 7,-4, 5,17,-1,  &
    -3, 8,-4, 6,18,-1,  &
    11,-5,13,-2,19, 1,  &
    12,10,14,-2,20, 2,  &
    -3,11,15,-2,21, 3,  &
    14,-5,16,10,22, 4,  &
    15,13,17,11,23, 5,  &
    -3,14,18,12,24, 6,  &
    17,-5,-4,13,25, 7,  &
    18,16,-4,14,26, 8,  &
    -3,17,-4,15,27, 9,  &
    20,-5,22,-2,-6,10,  &
    21,19,23,-2,-6,11,  &
    -3,20,24,-2,-6,12,  &
    23,-5,25,19,-6,13,  &
    24,22,26,20,-6,14,  &
    -3,23,27,21,-6,15,  &
    26,-5,-4,22,-6,16,  &
    27,25,-4,23,-6,17,  &
    -3,26,-4,24,-6,18   &
    ],pInt),[FE_NipNeighbors(me),FE_Nips(me)])


 ! *** FE_subNodeParent ***
 ! lists the group of nodes for which the center of gravity
 ! corresponds to the location of a each subnode.
 ! fill with 0.
 ! example: face-centered subnode with faceNodes 1,2,3,4 to be used in,
 !          e.g., a 8 IP grid, would be encoded:
 !          1, 2, 3, 4, 0, 0, 0, 0
 me = 0_pInt

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element   6 (2D 3node 1ip) has no subnodes
    0_pInt

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
    1, 2, 0,  & 
    2, 3, 0,  & 
    3, 1, 0,  &
    1, 2, 3   &
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])
 
 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
    1, 2, 0, 0,  & 
    2, 3, 0, 0,  & 
    3, 4, 0, 0,  &
    4, 1, 0, 0,  & 
    1, 2, 3, 4   &
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element  27 (2D 8node 9ip)
    reshape(int([&
    1, 1, 2, 0, 0, 0, 0, 0, 0,  &
    1, 2, 2, 0, 0, 0, 0, 0, 0,  &
    2, 2, 3, 0, 0, 0, 0, 0, 0,  &
    2, 3, 3, 0, 0, 0, 0, 0, 0,  &
    3, 3, 4, 0, 0, 0, 0, 0, 0,  &
    3, 4, 4, 0, 0, 0, 0, 0, 0,  &
    4, 4, 1, 0, 0, 0, 0, 0, 0,  &
    4, 1, 1, 0, 0, 0, 0, 0, 0,  &
    1, 1, 1, 1, 2, 2, 4, 4, 3,  &
    2, 2, 2, 2, 1, 1, 3, 3, 4,  &
    3, 3, 3, 3, 2, 2, 4, 4, 1,  &
    4, 4, 4, 4, 1, 1, 3, 3, 2   &
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element 134 (3D 4node 1ip) has no subnodes
    0_pInt

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
    1, 2, 0, 0,  & 
    2, 3, 0, 0,  & 
    3, 1, 0, 0,  &
    1, 4, 0, 0,  & 
    2, 4, 0, 0,  & 
    3, 4, 0, 0,  & 
    1, 2, 3, 0,  & 
    1, 2, 4, 0,  & 
    2, 3, 4, 0,  & 
    1, 3, 4, 0,  & 
    1, 2, 3, 4   &
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
    1, 2, 0, 0, 0, 0,  &
    2, 3, 0, 0, 0, 0,  &
    3, 1, 0, 0, 0, 0,  &
    1, 4, 0, 0, 0, 0,  &
    2, 5, 0, 0, 0, 0,  &
    3, 6, 0, 0, 0, 0,  &
    4, 5, 0, 0, 0, 0,  &
    5, 6, 0, 0, 0, 0,  &
    6, 4, 0, 0, 0, 0,  &
    1, 2, 3, 0, 0, 0,  &
    1, 2, 4, 5, 0, 0,  &
    2, 3, 5, 6, 0, 0,  &
    1, 3, 4, 6, 0, 0,  &
    4, 5, 6, 0, 0, 0,  &
    1, 2, 3, 4, 5, 6   &
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element 117 (3D 8node 1ip) has no subnodes
    0_pInt

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
    1, 2, 0, 0, 0, 0, 0, 0,  &
    2, 3, 0, 0, 0, 0, 0, 0,  &
    3, 4, 0, 0, 0, 0, 0, 0,  &
    4, 1, 0, 0, 0, 0, 0, 0,  &
    1, 5, 0, 0, 0, 0, 0, 0,  &
    2, 6, 0, 0, 0, 0, 0, 0,  &
    3, 7, 0, 0, 0, 0, 0, 0,  &
    4, 8, 0, 0, 0, 0, 0, 0,  &
    5, 6, 0, 0, 0, 0, 0, 0,  &
    6, 7, 0, 0, 0, 0, 0, 0,  &
    7, 8, 0, 0, 0, 0, 0, 0,  &
    8, 5, 0, 0, 0, 0, 0, 0,  &
    1, 2, 3, 4, 0, 0, 0, 0,  &
    1, 2, 6, 5, 0, 0, 0, 0,  &
    2, 3, 7, 6, 0, 0, 0, 0,  &
    3, 4, 8, 7, 0, 0, 0, 0,  &
    1, 4, 8, 5, 0, 0, 0, 0,  &
    5, 6, 7, 8, 0, 0, 0, 0,  &
    1, 2, 3, 4, 5, 6, 7, 8   & 
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])

 me = me + 1_pInt
 FE_subNodeParent(1:FE_Nips(me),1:FE_NsubNodes(me),me) = &  ! element  21 (3D 20node 27ip)
    reshape(int([&
    1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    1, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 4, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    1, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    5, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    5, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 2, 2, 4, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 2, 2, 2, 1, 1, 3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    3, 3, 3, 3, 2, 2, 4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    4, 4, 4, 4, 1, 1, 3, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 2, 2, 5, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 2, 2, 2, 1, 1, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 2, 2, 2, 3, 3, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    3, 3, 3, 3, 2, 2, 7, 7, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    3, 3, 3, 3, 4, 4, 7, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    4, 4, 4, 4, 3, 3, 8, 8, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    4, 4, 4, 4, 1, 1, 8, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 4, 4, 5, 5, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    5, 5, 5, 5, 1, 1, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 6, 6, 6, 2, 2, 5, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 6, 6, 6, 2, 2, 7, 7, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 7, 7, 3, 3, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 7, 7, 3, 3, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 8, 8, 4, 4, 7, 7, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 8, 8, 4, 4, 5, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    5, 5, 5, 5, 1, 1, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    5, 5, 5, 5, 6, 6, 8, 8, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 6, 6, 6, 5, 5, 7, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 7, 7, 6, 6, 8, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 8, 8, 5, 5, 7, 7, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 5, 5, 5, 5, 3, 3, 6, 6, 8, 8, 7, &
    2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 3, 6, 6, 6, 6, 4, 4, 5, 5, 7, 7, 8, &
    3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4, 7, 7, 7, 7, 1, 1, 6, 6, 8, 8, 5, &
    4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 3, 3, 3, 3, 8, 8, 8, 8, 2, 2, 5, 5, 7, 7, 6, &
    5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 6, 6, 6, 6, 8, 8, 8, 8, 2, 2, 4, 4, 7, 7, 3, &
    6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 5, 5, 5, 5, 7, 7, 7, 7, 1, 1, 3, 3, 8, 8, 4, &
    7, 7, 7, 7, 7, 7, 7, 7, 3, 3, 3, 3, 6, 6, 6, 6, 8, 8, 8, 8, 2, 2, 4, 4, 5, 5, 1, &
    8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4, 4, 5, 5, 5, 5, 7, 7, 7, 7, 1, 1, 3, 3, 6, 6, 2  & 
    ],pInt),[FE_Nips(me),FE_NsubNodes(me)])


 ! *** FE_subNodeOnIPFace ***
 ! indicates which subnodes make up the interfaces enclosing the IP volume.
 ! The sorting convention is such that the outward pointing normal
 ! follows from a right-handed traversal of the face node list.
 me = 0_pInt
 
 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
     2, 3, & ! 1
     3, 1, &
     1, 2  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
     4, 7, & ! 1
     6, 1, &
     7, 6, &
     1, 4, &
     2, 5, & ! 2
     7, 4, &
     5, 7, &
     4, 2, &
     7, 5, & ! 3
     3, 6, &
     5, 3, &
     6, 7  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
     5, 9, & ! 1
     8, 1, &
     9, 8, &
     1, 5, &
     2, 6, & ! 2
     9, 5, &
     6, 9, &
     5, 2, &
     9, 7, & ! 3
     4, 8, &
     7, 4, &
     8, 9, &
     6, 3, & ! 4
     7, 9, &
     3, 7, &
     9, 6  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element  27 (2D 8node 9ip)
    reshape(int([&
     5,13, & ! 1
    12, 1, &
    13,12, &
     1, 5, &
     6,14, & ! 2
    13, 5, &
    14,13, &
     5, 6, &
     2, 7, & ! 3
    14, 6, &
     7,14, &
     6, 2, &
    13,16, & ! 4
    11,12, &
    16,11, &
    12,13, &
    14,15, & ! 5
    16,13, &
    15,16, &
    13,14, &
     7, 8, & ! 6
    15,14, &
     8,15, &
    14, 7, &
    16,10, & ! 7
     4,11, &
    10, 4, &
    11,16, &
    15, 9, & ! 8
    10,16, &
     9,10, &
    16,15, &
     8, 3, & ! 9
     9,15, &
     3, 9, &
    15, 8  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
     1, 3, 2, & ! 1
     1, 2, 4, &
     2, 3, 4, &
     1, 4, 3  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
     5,11,15,12 , & ! 1
     1, 8,14, 7 , &
     7,14,15,11 , &
     1, 5,12, 8 , &
     8,12,15,14 , &
     1, 7,11, 5 , &
     2, 6,13, 9 , & ! 2
     5,12,15,11 , &
     6,11,15,13 , &
     2, 9,12, 5 , &
     9,13,15,12 , &
     2,13,11, 6 , &
     6,13,15,11 , & ! 3
     3, 7,14,10 , &
     3,10,13, 6 , &
     7,11,15,14 , &
    13,10,14,15 , &
     3, 6,11, 7 , &
     9,12,15,13 , & ! 4
     4,10,14, 8 , &
    10,13,15,14 , &
     4, 8,12, 9 , &
     4, 9,13,10 , &
     8,10,15,12   &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
     7,16,21,17, & ! 1
     1,10,19, 9, &
     9,19,21,16, &
     1, 7,17,10, &
    10,17,21,19, &
     1, 9,16, 7, &
     2, 8,18,11, & ! 2
     7,17,21,16, &
     8,16,21,18, &
     2,11,17, 7, &
    11,18,21,17, &
     2, 7,16, 8, &
     8,18,21,16, & ! 3
     3, 9,19,12, &
     3,12,18, 8, &
     9,16,21,19, &
    12,19,21,18, &
     3, 8,16, 9, &
    13,17,21,20, & ! 4
     4,15,19,10, &
    15,20,21,19, &
     4,10,17,13, &
     4,13,20,15, &
    10,19,21,17, &
     5,11,18,14, & ! 5
    13,20,21,17, &
    14,18,21,20, &
     5,13,17,11, &
     5,14,20,13, &
    11,17,21,18, &
    14,20,21,18, & ! 6
     6,12,19,15, &
     6,14,18,12, &
    15,19,21,20, &
     6,15,20,14, &
    12,18,21,19  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
     2, 3, 7, 6, & ! 1
     1, 5, 8, 4, &
     3, 4, 8, 7, &
     1, 2, 6, 5, &
     5, 6, 7, 8, &
     1, 4, 3, 2  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
     9,21,27,22, & ! 1
     1,13,25,12, &
    12,25,27,21, &
     1, 9,22,13, &
    13,22,27,25, &
     1,12,21, 9, &
     2,10,23,14, & ! 2
     9,22,27,21, &
    10,21,27,23, &
     2,14,22, 9, &
    14,23,27,22, &
     2, 9,21,10, &
    11,24,27,21, & ! 3
     4,12,25,16, &
     4,16,24,11, &
    12,21,27,25, &
    16,25,27,24, &
     4,11,21,12, &
     3,15,23,10, & ! 4
    11,21,27,24, &
     3,11,24,15, &
    10,23,27,21, &
    15,24,27,23, &
     3,10,21,11, &
    17,22,27,26, & ! 5
     5,20,25,13, &
    20,26,27,25, &
     5,13,22,17, &
     5,17,26,20, &
    13,25,27,22, &
     6,14,23,18, & ! 6
    17,26,27,22, &
    18,23,27,26, &
     6,17,22,14, &
     6,18,26,17, &
    14,22,27,23, &
    19,26,27,24, & ! 7
     8,16,25,20, &
     8,19,24,16, &
    20,25,27,26, &
     8,20,26,19, &
    16,24,27,25, &
     7,18,23,15, & ! 8
    19,24,27,26, &
     7,15,24,19, &
    18,26,27,23, &
     7,19,26,18, &
    15,23,27,24  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_subNodeOnIPFace(1:FE_NipFaceNodes(me),1:FE_NipNeighbors(me),1:FE_Nips(me),me) = &  ! element  21 (3D 20node 27ip)
    reshape(int([&
     9,33,57,37, &  ! 1
     1,17,44,16, &
    33,16,44,57, &
     1, 9,37,17, &
    17,37,57,44, &
     1,16,33, 9, &
    10,34,58,38, &  !  2
     9,37,57,33, &
    34,33,57,58, &
     9,10,38,37, &
    37,38,58,57, &
     9,33,34,10, &
     2,11,39,18, &  !  3
    10,38,58,34, &
    11,34,58,39, &
    10, 2,18,38, &
    38,18,39,58, &
    10,34,11, 2, &
    33,36,60,57, &  !  4
    16,44,43,15, &
    36,15,43,60, &
    16,33,57,44, &
    44,57,60,43, &
    16,15,36,33, &
    34,35,59,58, &  !  5
    33,57,60,36, &
    35,36,60,59, &
    33,34,58,57, &
    57,58,59,60, &
    33,36,35,34, &
    11,12,40,39, &  !  6
    34,58,59,35, &
    12,35,59,40, &
    34,11,39,58, &
    58,39,40,59, &
    34,35,12,11, &
    36,14,42,60, &  !  7
    15,43,20, 4, &
    14, 4,20,42, &
    15,36,60,43, &
    43,60,42,20, &
    15, 4,14,36, &
    35,13,41,59, &  !  8
    36,60,42,14, &
    13,14,42,41, &
    36,35,59,60, &
    60,59,41,42, &
    36,14,13,35, &
    12, 3,19,40, &  !  9
    35,59,41,13, &
     3,13,41,19, &
    35,12,40,59, &
    59,40,19,41, &
    35,13, 3,12, &
    37,57,61,45, &  ! 10
    17,21,52,44, &
    57,44,52,61, &
    17,37,45,21, &
    21,45,61,52, &
    17,44,57,37, &
    38,58,62,46, &  ! 11
    37,45,61,57, &
    58,57,61,62, &
    37,38,46,45, &
    45,46,62,61, &
    37,57,58,38, &
    18,39,47,22, &  ! 12
    38,46,62,58, &
    39,58,62,47, &
    38,18,22,46, &
    46,22,47,62, &
    38,58,39,18, &
    57,60,64,61, &  ! 13
    44,52,51,43, &
    60,43,51,64, &
    44,57,61,52, &
    52,61,64,51, &
    44,43,60,57, &
    58,59,63,62, &  ! 14
    57,61,64,60, &
    59,60,64,63, &
    57,58,62,61, &
    61,62,63,64, &
    57,60,59,58, &
    39,40,48,47, &  ! 15
    58,62,63,59, &
    40,59,63,48, &
    58,39,47,62, &
    62,47,48,63, &
    58,59,40,39, &
    60,42,50,64, &  ! 16
    43,51,24,20, &
    42,20,24,50, &
    43,60,64,51, &
    51,64,50,24, &
    43,20,42,60, &
    59,41,49,63, &  ! 17
    60,64,50,42, &
    41,42,50,49, &
    60,59,63,64, &
    64,63,49,50, &
    60,42,41,59, &
    40,19,23,48, &  ! 18
    59,63,49,41, &
    19,41,49,23, &
    59,40,48,63, &
    63,48,23,49, &
    59,41,19,40, &
    45,61,53,25, &  ! 19
    21, 5,32,52, &
    61,52,32,53, &
    21,45,25, 5, &
     5,25,53,32, &
    21,52,61,45, &
    46,62,54,26, &  ! 20
    45,25,53,61, &
    62,61,53,54, &
    45,46,26,25, &
    25,26,54,53, &
    45,61,62,46, &
    22,47,27, 6, &  ! 21
    46,26,54,62, &
    47,62,54,27, &
    46,22, 6,26, &
    26, 6,27,54, &
    46,62,47,22, &
    61,64,56,53, &  ! 22
    52,32,31,51, &
    64,51,31,56, &
    52,61,53,32, &
    32,53,56,31, &
    52,51,64,61, &
    62,63,55,54, &  ! 23
    61,53,56,64, &
    63,64,56,55, &
    61,62,54,53, &
    53,54,55,56, &
    61,64,63,62, &
    47,48,28,27, &  ! 24
    62,54,55,63, &
    48,63,55,28, &
    62,47,27,54, &
    54,27,28,55, &
    62,63,48,47, &
    64,50,30,56, &  ! 25
    51,31, 8,24, &
    50,24, 8,30, &
    51,64,56,31, &
    31,56,30, 8, &
    51,24,50,64, &
    63,49,29,55, &  ! 26
    64,56,30,50, &
    49,50,30,29, &
    64,63,55,56, &
    56,55,29,30, &
    64,50,49,63, &
    48,23, 7,28, &  ! 27
    63,55,29,49, &
    23,49,29, 7, &
    63,48,28,55, &
    55,28, 7,29, &
    63,49,23,48  &
    ],pInt),[FE_NipFaceNodes(me),FE_NipNeighbors(me),FE_Nips(me)])


end subroutine mesh_build_FEdata

end module mesh

