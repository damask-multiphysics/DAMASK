!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant damage field
!--------------------------------------------------------------------------------------------------
module damage_none
  use config
  use material

  implicit none
  public

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine damage_none_init

  integer :: h,Nmaterialpoints

  print'(/,a)', ' <<<+-  damage_none init  -+>>>'; flush(6)

  do h = 1, size(material_name_homogenization)
    if (damage_type(h) /= DAMAGE_NONE_ID) cycle

    Nmaterialpoints = count(material_homogenizationAt == h)
    damageState(h)%sizeState = 0
    allocate(damageState(h)%state0   (0,Nmaterialpoints))
    allocate(damageState(h)%subState0(0,Nmaterialpoints))
    allocate(damageState(h)%state    (0,Nmaterialpoints))

    deallocate(damage(h)%p)
    allocate  (damage(h)%p(1), source=damage_initialPhi(h))

  enddo

end subroutine damage_none_init

end module damage_none
