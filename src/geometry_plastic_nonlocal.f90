!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Geometric information about the IP cells needed for the nonlocal
! plasticity model
!--------------------------------------------------------------------------------------------------
module geometry_plastic_nonlocal
  use prec
  use results

  implicit none
  public

  integer, protected :: &
    geometry_plastic_nonlocal_nIPneighbors

  integer,     dimension(:,:,:,:), allocatable, protected :: &
    geometry_plastic_nonlocal_IPneighborhood                                                        !< 6 or less neighboring IPs as [element ID, IP ID, face ID that point to me]

  real(pReal), dimension(:,:),     allocatable, protected :: &
    geometry_plastic_nonlocal_IPvolume0                                                             !< volume associated with IP (initially!)

  real(pReal), dimension(:,:,:),   allocatable, protected :: &
    geometry_plastic_nonlocal_IParea0                                                               !< area of interface to neighboring IP (initially!)

  real(pReal), dimension(:,:,:,:), allocatable, protected :: &
    geometry_plastic_nonlocal_IPareaNormal0                                                         !< area normal of interface to neighboring IP (initially!)


contains

!---------------------------------------------------------------------------------------------------
!> @brief Set the integration point (IP) neighborhood
!> @details: The IP neighborhood for element ID (last index), IP ID (second but last index) and
!            face ID (second index) gives the element ID (1 @ first index), IP ID (2 @ first index)
!            and face ID (3 @ first index).
!            A triangle (2D) has 3 faces, a quadrilateral (2D) had 4 faces, a tetrahedron (3D) has
!            4 faces, and a hexahedron (3D) has 6 faces.
!---------------------------------------------------------------------------------------------------
subroutine geometry_plastic_nonlocal_setIPneighborhood(IPneighborhood)

  integer, dimension(:,:,:,:), intent(in) :: IPneighborhood

  geometry_plastic_nonlocal_IPneighborhood = IPneighborhood
  geometry_plastic_nonlocal_nIPneighbors   = size(IPneighborhood,2)


end subroutine geometry_plastic_nonlocal_setIPneighborhood


!---------------------------------------------------------------------------------------------------
!> @brief Set the initial volume associated with an integration point
!---------------------------------------------------------------------------------------------------
subroutine geometry_plastic_nonlocal_setIPvolume(IPvolume)

  real(pReal), dimension(:,:), intent(in) :: IPvolume

  geometry_plastic_nonlocal_IPvolume0 = IPvolume

end subroutine geometry_plastic_nonlocal_setIPvolume


!---------------------------------------------------------------------------------------------------
!> @brief Set the initial areas of the unit triangle/unit quadrilateral/tetrahedron/hexahedron
!         encompassing an integration point
!---------------------------------------------------------------------------------------------------
subroutine geometry_plastic_nonlocal_setIParea(IParea)

  real(pReal), dimension(:,:,:), intent(in) :: IParea

  geometry_plastic_nonlocal_IParea0 = IParea

end subroutine geometry_plastic_nonlocal_setIParea


!---------------------------------------------------------------------------------------------------
!> @brief Set the direction normal of the areas of the triangle/quadrilateral/tetrahedron/hexahedron
!         encompassing an integration point
!---------------------------------------------------------------------------------------------------
subroutine geometry_plastic_nonlocal_setIPareaNormal(IPareaNormal)

  real(pReal), dimension(:,:,:,:), intent(in) :: IPareaNormal

  geometry_plastic_nonlocal_IPareaNormal0 = IPareaNormal

end subroutine geometry_plastic_nonlocal_setIPareaNormal


!---------------------------------------------------------------------------------------------------
!> @brief Free memory used by variables only needed by plastic_nonlocal
!---------------------------------------------------------------------------------------------------
subroutine geometry_plastic_nonlocal_disable

  if(allocated(geometry_plastic_nonlocal_IPneighborhood)) &
    deallocate(geometry_plastic_nonlocal_IPneighborhood)

  if(allocated(geometry_plastic_nonlocal_IPvolume0)) &
    deallocate(geometry_plastic_nonlocal_IPvolume0)

  if(allocated(geometry_plastic_nonlocal_IParea0)) &
    deallocate(geometry_plastic_nonlocal_IParea0)

  if(allocated(geometry_plastic_nonlocal_IPareaNormal0)) &
    deallocate(geometry_plastic_nonlocal_IPareaNormal0)

end subroutine geometry_plastic_nonlocal_disable


!---------------------------------------------------------------------------------------------------
!> @brief Write geometry data to results file
!---------------------------------------------------------------------------------------------------
subroutine geometry_plastic_nonlocal_results

  integer,     dimension(:),   allocatable :: shp

  call results_openJobFile

  writeVolume: block
    real(pReal), dimension(:), allocatable :: temp
    shp = shape(geometry_plastic_nonlocal_IPvolume0)
    temp = reshape(geometry_plastic_nonlocal_IPvolume0,[shp(1)*shp(2)])
    call results_writeDataset(temp,'geometry','v_0',&
                              'initial cell volume','m³')
  end block writeVolume

  writeAreas: block
    real(pReal), dimension(:,:), allocatable :: temp
    shp = shape(geometry_plastic_nonlocal_IParea0)
    temp = reshape(geometry_plastic_nonlocal_IParea0,[shp(1),shp(2)*shp(3)])
    call results_writeDataset(temp,'geometry','a_0',&
                              'initial cell face area','m²')
  end block writeAreas

  writeNormals: block
    real(pReal), dimension(:,:,:), allocatable :: temp
    shp = shape(geometry_plastic_nonlocal_IPareaNormal0)
    temp = reshape(geometry_plastic_nonlocal_IPareaNormal0,[shp(1),shp(2),shp(3)*shp(4)])
    call results_writeDataset(temp,'geometry','n_0',&
                              'initial cell face normals','-',transposed=.false.)
  end block writeNormals


  call results_closeJobFile

end subroutine geometry_plastic_nonlocal_results

end module geometry_plastic_nonlocal
