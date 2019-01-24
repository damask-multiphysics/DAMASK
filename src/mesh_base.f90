!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc,FEM, Abaqus and the spectral solver
!--------------------------------------------------------------------------------------------------
module mesh_base

 use, intrinsic :: iso_c_binding
 use prec, only: &
   pStringLen, &
   pReal, &
   pInt
 use element, only: &
   tElement

  implicit none

!---------------------------------------------------------------------------------------------------
!> Properties of a the whole mesh (consisting of one type of elements)
!---------------------------------------------------------------------------------------------------
 type, public :: tMesh
   type(tElement) :: &
     elem
   real(pReal), dimension(:,:), allocatable, public :: &
     ipVolume, &                                                                                 !< volume associated with each IP (initially!)
     node0, &                                                                                    !< node x,y,z coordinates (initially)
     node                                                                                        !< node x,y,z coordinates (deformed)
   integer(pInt), dimension(:,:), allocatable, public :: &    
     cellnodeParent                                                                               !< cellnode's parent element ID, cellnode's intra-element ID 
   character(pStringLen) :: solver = "undefined"
   integer(pInt)         :: &
     Nnodes, &                                                                                   !< total number of nodes in mesh
     Nelems = -1_pInt, &
     elemType, &
     Ncells, &
     nIPneighbors, &
     NcellNodes, &
     maxElemsPerNode
   integer(pInt), dimension(:), allocatable, public :: &
     homogenizationAt, &
     microstructureAt
   integer(pInt), dimension(:,:), allocatable, public :: &
     connectivity
 end type tMesh

end module mesh_base
