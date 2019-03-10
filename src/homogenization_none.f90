!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief dummy homogenization homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
module homogenization_none

  implicit none
  private
  
  public :: &
    homogenization_none_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_none_init()
  use debug, only: &
    debug_HOMOGENIZATION, &
    debug_level, &
    debug_levelBasic
  use config, only: &
    config_homogenization
  use material, only: &
    homogenization_type, &
    material_homogenizationAt, &
    homogState, &
    HOMOGENIZATION_NONE_LABEL, &
    HOMOGENIZATION_NONE_ID
 
  implicit none
  integer :: &
    Ninstance, &
    h, &
    NofMyHomog
 
  write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_NONE_label//' init  -+>>>'
 
  Ninstance = count(homogenization_type == HOMOGENIZATION_NONE_ID)
  if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
  do h = 1, size(homogenization_type)
    if (homogenization_type(h) /= HOMOGENIZATION_NONE_ID) cycle
    
    NofMyHomog = count(material_homogenizationAt == h)
    homogState(h)%sizeState = 0
    homogState(h)%sizePostResults = 0
    allocate(homogState(h)%state0   (0,NofMyHomog))
    allocate(homogState(h)%subState0(0,NofMyHomog))
    allocate(homogState(h)%state    (0,NofMyHomog))
 
  enddo

end subroutine homogenization_none_init

end module homogenization_none
