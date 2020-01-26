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

  integer :: h,NofMyHomog

  write(6,'(/,a)') ' <<<+-  damage_'//DAMAGE_NONE_LABEL//' init  -+>>>'; flush(6)

  do h = 1, size(config_homogenization)
    if (damage_type(h) /= DAMAGE_NONE_ID) cycle

    NofMyHomog = count(material_homogenizationAt == h)
    damageState(h)%sizeState = 0
    allocate(damageState(h)%state0   (0,NofMyHomog))
    allocate(damageState(h)%subState0(0,NofMyHomog))
    allocate(damageState(h)%state    (0,NofMyHomog))
    
    deallocate(damage(h)%p)
    allocate  (damage(h)%p(1), source=damage_initialPhi(h))
      
  enddo

end subroutine damage_none_init

end module damage_none
