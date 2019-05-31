!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant damage field
!--------------------------------------------------------------------------------------------------
module damage_none
  use config
  use material

  implicit none
  private
  
  public :: &
    damage_none_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine damage_none_init

  integer :: &
    homog, &
    NofMyHomog

  write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_NONE_LABEL//' init  -+>>>'

   initializeInstances: do homog = 1, size(config_homogenization)
    
    myhomog: if (damage_type(homog) == DAMAGE_NONE_ID) then
      NofMyHomog = count(material_homogenizationAt == homog)
      damageState(homog)%sizeState = 0
      damageState(homog)%sizePostResults = 0
      allocate(damageState(homog)%state0   (0,NofMyHomog))
      allocate(damageState(homog)%subState0(0,NofMyHomog))
      allocate(damageState(homog)%state    (0,NofMyHomog))
      
      deallocate(damage(homog)%p)
      allocate  (damage(homog)%p(1), source=damage_initialPhi(homog))
      
    endif myhomog
  enddo initializeInstances

end subroutine damage_none_init

end module damage_none
